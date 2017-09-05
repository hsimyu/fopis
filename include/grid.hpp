#ifndef __TDPIC_GRID_H_INCLUDED__
#define __TDPIC_GRID_H_INCLUDED__
#include <vector>
#include <map>
#include "particle.hpp"
#include "environment.hpp"
#include "spacecraft.hpp"

#define H5_USE_BOOST
#include <highfive/H5File.hpp>

class ChildGrid;

//! @class Grid
class Grid  : public std::enable_shared_from_this<Grid> {
    protected:
        //! グリッド関係ツリー
        std::vector<std::unique_ptr<ChildGrid>> children;

        // オブジェクト定義は Node ベース
        std::vector<Spacecraft> objects;

        //! Particleを格納したstd::vectorを、粒子種ごとに保持したstd::vector
        ParticleArray particles;

        //! 親のどの座標にくっついているか
        //! @{
        int from_ix;
        int from_iy;
        int from_iz;
        int to_ix;
        int to_iy;
        int to_iz;
        double base_x;
        double base_y;
        double base_z;
        //! @}

        unsigned int id;
        int nx, ny, nz;
        int level;
        int sumTotalNumOfChildGrids;
        double dx;
        double dt;

        std::unique_ptr<Field> field;

        // Class Unique ID
        static unsigned int nextID;

        // データ出力用のnodes/edges/faces/cellsを取り出す
        boost::multi_array<float, 3> getTrueNodes(const tdArray& x3D, const double unnorm = 1.0) const;
        boost::multi_array<float, 3> getTrueNodes(const RhoArray& rho, const int pid, const double unnorm = 1.0) const;
        boost::multi_array<float, 3> getDensity(const int) const;

        //! データ出力用のnodesの数を返す
        virtual int getXNodeSize(void) const = 0;
        virtual int getYNodeSize(void) const = 0;
        virtual int getZNodeSize(void) const = 0;

        // 粒子更新用のメソッド
        //! ES: 静電
        //! EM: 電磁
        void updateParticleVelocityES(void);
        void updateParticlePositionES(void);
        void updateParticleVelocityEM(void);
        void updateParticlePositionEM(void);

        //! 場の初期化
        void initializeField();

        //! ポアソンソルバの実体
        virtual void solvePoissonPSOR(const int loopnum) = 0;
        virtual double checkPhiResidual() = 0;
    public:
        Grid();
        ~Grid();

        unsigned int getNextID(void) {
            return nextID++;
        }

        // IDへのsetterは提供しない、代わりにresetterを用意
        void resetNextID(void) {
            nextID = 0;
        }
        unsigned int getID(void) const { return id; }

        // accessors
        void   setFromIX(int _ix) { from_ix = _ix; }
        void   setFromIY(int _iy) { from_iy = _iy; }
        void   setFromIZ(int _iz) { from_iz = _iz; }
        int    getFromIX(void) const { return from_ix; }
        int    getFromIY(void) const { return from_iy; }
        int    getFromIZ(void) const { return from_iz; }

        void   setToIX(int _ix) { to_ix = _ix; }
        void   setToIY(int _iy) { to_iy = _iy; }
        void   setToIZ(int _iz) { to_iz = _iz; }
        int    getToIX(void) const { return to_ix; }
        int    getToIY(void) const { return to_iy; }
        int    getToIZ(void) const { return to_iz; }

        void   setBaseX(double _x){ base_x = _x; }
        void   setBaseY(double _y){ base_y = _y; }
        void   setBaseZ(double _z){ base_z = _z; }
        double getBaseX(void)  const { return base_x; }
        double getBaseY(void)  const { return base_y; }
        double getBaseZ(void)  const { return base_z; }

        void setNX(int _x){ nx = _x; }
        void setNY(int _y){ ny = _y; }
        void setNZ(int _z){ nz = _z; }
        int  getNX(void) const { return nx; }
        int  getNY(void) const { return ny; }
        int  getNZ(void) const { return nz; }

        void setLevel(int l){ level = l; }
        int  getLevel(void) const { return level; }
        void   setDx(double _dx){ dx = _dx; }
        double getDx(void) const { return dx; }
        void   setDt(double _dt){ dt = _dt; }
        double getDt(void) const { return dt; }

        //! 場
        tdArray& getPhi() { return field->getPhi(); }
        RhoArray& getRho() { return field->getRho(); }

        //! 物体関連
        const std::vector<Spacecraft>& getObjects() const { return objects; };

        //! 粒子関連
        ParticleArray& getParticles() {return particles;}
        size_t getValidParticleNumber(const int) const;

        // 親子でのScalarやりとり用
        void copyScalarToChildren(std::string);
        void updateChildrenPhi();
        void correctChildrenPhi();
        void interpolateRhoValueToChildren();
        void restrictPhiValueToChildren();

        // 子供管理メソッド
        void makeChild(const int, const int, const int, const int, const int, const int);
        void addChild(std::unique_ptr<ChildGrid>&&);
        void removeChild(const int);

        int getChildrenLength() const {
            return children.size();
        }

        // child gridの総数を保存するためのメソッド
        int getSumOfChild(void) const {
            return sumTotalNumOfChildGrids;
        }

        void setSumOfChild(const int s) {
            sumTotalNumOfChildGrids = s;
        }

        virtual void incrementSumOfChild(void) = 0;
        virtual void decrementSumOfChild(void) = 0;

        // update fields は各Gridクラスで実装する
        virtual void updateRho(void) = 0;
        virtual void solvePoisson(void) = 0;
        virtual void updateEfield(void) = 0;
        virtual void updateEfieldFDTD(void) = 0;
        virtual void updateBfield(void) = 0;

        // update particles
        void updateParticleVelocity(void);
        void updateParticlePosition(void);

        //! エネルギーを計算する
        double getParticleEnergy(void) const;
        double getEFieldEnergy(void) const;
        double getBFieldEnergy(void) const;

        // create mesh nodes array
        float** getMeshNodes(int);
        int* getChildOfPatches(void);
        int* getNumOfPatches(void);
        int* getChildIdMap(void);
        int getMaxLevel(void);

        int getMaxChildLevel() {
            return this->getMaxLevel() - level;
        }

        // {id: [子のIDを格納したvector]}のmapを作成する
        std::map<int, std::vector<int> > getChildMapOnRoot(void);
        void addChildrenIDToMap(std::map<int, std::vector<int> >&);
        std::vector<int> getChildrenIDs(void);

        // [level: [IDを格納したvector]]のvectorを作成する
        std::vector< std::vector<int> > getIDMapOnRoot(void);
        void addIDToVector(std::vector< std::vector<int> >&);

        // HDF5にデータを突っ込む
        void putFieldData(HighFive::Group& group, const std::string& data_type_name, const std::string& i_timestamp) const;

        void printInfo() const;

        //! 粒子境界チェック
        virtual void checkXBoundary(ParticleArray& pbuff, Particle& p, const double slx) = 0;
        virtual void checkYBoundary(ParticleArray& pbuff, Particle& p, const double sly) = 0;
        virtual void checkZBoundary(ParticleArray& pbuff, Particle& p, const double slz) = 0;
};

class RootGrid : public Grid {
    public:
        RootGrid();

        void injectParticlesFromBoundary(void);

        //! 物体定義初期化
        void initializeObject();
        void resetObjects();

        //! 物体のキャパシタンス行列初期化または読み込み
        void initializeObjectsCmatrix(void);
        void emitParticlesFromObjects(void);

        //! 場の初期化 / 更新
        virtual void updateRho(void) override;
        virtual void solvePoisson(void) override;
        virtual void updateEfield(void) override;
        virtual void updateEfieldFDTD(void) override;
        virtual void updateBfield(void) override;

        virtual void incrementSumOfChild(void) override;
        virtual void decrementSumOfChild(void) override;

        virtual void checkXBoundary(ParticleArray& pbuff, Particle& p, const double slx) override {
            if(p.x < 0.0) {
                if (Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low)) {
                    // 計算空間の端でない場合は粒子を隣へ送る
                    // 計算空間の端にいるが、周期境界の場合も粒子を送る必要がある -> isNotBoundary()でまとめて判定できる
                    pbuff[0].push_back(p);
                }
                p.makeInvalid();
            } else if (p.x > (slx - dx)) {
                if (Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up)) {
                    //! 計算空間の端でない場合、slx - dx から slx までの空間は下側の空間が担当するため、 slx を超えた場合のみ粒子を送信する
                    //! また、計算空間の端にいるが、周期境界の場合も同様の処理でよいため、isNotBoundary()でまとめて判定できる
                    if (p.x > slx) {
                        pbuff[1].push_back(p);
                        p.makeInvalid();
                    }
                } else {
                    p.makeInvalid();
                }
            }
        }

        virtual void checkYBoundary(ParticleArray& pbuff, Particle& p, const double sly) override {
            if(p.y < 0.0) {
                if (Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low)) pbuff[2].push_back(p);
                p.makeInvalid();
            } else if (p.y > (sly - dx)) {
                if (Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up)) {
                    if (p.y > sly) {
                        pbuff[3].push_back(p);
                        p.makeInvalid();
                    }
                } else {
                    p.makeInvalid();
                }
            }
        }

        virtual void checkZBoundary(ParticleArray& pbuff, Particle& p, const double slz) override {
            if(p.z < 0.0) {
                if (Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low)) pbuff[4].push_back(p);
                p.makeInvalid();
            } else if (p.z > (slz - dx)) {
                if (Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up)) {
                    if (p.z > slz) {
                        pbuff[5].push_back(p);
                        p.makeInvalid();
                    }
                } else {
                    p.makeInvalid();
                }
            }
        }

    private:
        //! 基本的には root_grid 中に対象の点(Object定義点)が含まれているかを判定するために呼ぶ
        //! i, j, k は整数座標(全体の計算空間上の)
        bool isInnerNode(const int i, const int j, const int k) const {
            return 
                (i >= Environment::getAssignedXBegin()) &&
                (i <= Environment::getAssignedXEnd()) &&
                (j >= Environment::getAssignedYBegin()) &&
                (j <= Environment::getAssignedYEnd()) &&
                (k >= Environment::getAssignedZBegin()) &&
                (k <= Environment::getAssignedZEnd());
        }

        //! Glueノード込みで内部かどうかを判定する
        bool isInnerNodeWithGlue(const int i, const int j, const int k) const {
            return 
                (i >= (Environment::getAssignedXBegin() - 1)) &&
                (i <= (Environment::getAssignedXEnd() + 1)) &&
                (j >= (Environment::getAssignedYBegin() - 1)) &&
                (j <= (Environment::getAssignedYEnd() + 1)) &&
                (k >= (Environment::getAssignedZBegin() - 1)) &&
                (k <= (Environment::getAssignedZEnd() + 1));
        }

        //! Glueノードであるかどうかを判定する
        bool isGlueNode(const int i, const int j, const int k) const {
            return isInnerNodeWithGlue(i, j, k) && !isInnerNode(i, j, k);
        }

        //! Glueセル込みで内部かどうかを判定する
        //! セル周りのノードが全て内部(Glue含)ならInnerと判断できる
        bool isInnerCellWithGlue(const int i, const int j, const int k) const {
            return
                isInnerNodeWithGlue(i,   j,   k) &&
                isInnerNodeWithGlue(i+1, j,   k) &&
                isInnerNodeWithGlue(i,   j+1, k) &&
                isInnerNodeWithGlue(i+1, j+1, k) &&
                isInnerNodeWithGlue(i,   j,   k+1) &&
                isInnerNodeWithGlue(i+1, j,   k+1) &&
                isInnerNodeWithGlue(i,   j+1, k+1) &&
                isInnerNodeWithGlue(i+1, j+1, k+1);
        }

        virtual int getXNodeSize(void) const override;
        virtual int getYNodeSize(void) const override;
        virtual int getZNodeSize(void) const override;

        virtual void solvePoissonPSOR(const int loopnum) override;
        virtual double checkPhiResidual() override;
};

class ChildGrid : public Grid {
    public:
        ChildGrid(std::shared_ptr<Grid>, const int, const int, const int, const int, const int, const int);

        virtual void updateRho(void) override;
        virtual void solvePoisson(void) override;
        void solvePoissonPost(void);
        virtual void updateEfield(void) override;
        virtual void updateEfieldFDTD(void) override;
        virtual void updateBfield(void) override;

        void  setParent(std::shared_ptr<Grid> g){ parent = g; }
        std::shared_ptr<Grid> getParent(){ return parent; }
        void copyPhiToParent();
        virtual void incrementSumOfChild(void) override;
        virtual void decrementSumOfChild(void) override;

    private:
        std::shared_ptr<Grid> parent;
        void checkGridValidness();
        virtual int getXNodeSize(void) const override;
        virtual int getYNodeSize(void) const override;
        virtual int getZNodeSize(void) const override;
        virtual void checkXBoundary(ParticleArray& pbuff, Particle& p, const double slx) override;
        virtual void checkYBoundary(ParticleArray& pbuff, Particle& p, const double sly) override;
        virtual void checkZBoundary(ParticleArray& pbuff, Particle& p, const double slz) override;
        virtual void solvePoissonPSOR(const int loopnum) override;
        virtual double checkPhiResidual() override;
};
#endif
