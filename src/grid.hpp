#ifndef __TDPIC_GRID_H_INCLUDED__
#define __TDPIC_GRID_H_INCLUDED__
#include <vector>
#include <map>
#include <silo.h>
#include "particle.hpp"
#include "field.hpp"
#include "environment.hpp"
#include "utils.hpp"
#include "spacecraft.hpp"

//! @class Grid
class Grid {
    private:
        //! グリッド関係ツリー
        Grid* parent;
        std::vector<Grid*> children;

        // オブジェクト定義は Node ベース
        std::vector<Spacecraft> objects;

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
        float* getTrueNodes(const tdArray&);

        //! データ出力用のnodesの数を返す
        int getXNodeSize(void) const;
        int getYNodeSize(void) const;
        int getZNodeSize(void) const;

        // 粒子更新用のメソッド
        //! ES: 静電
        //! EM: 電磁
        void updateParticleVelocityES(void);
        void updateParticlePositionES(void);
        void updateParticleVelocityEM(void);
        void updateParticlePositionEM(void);

    public:
        Grid(void);
        Grid(Grid*, const int, const int, const int, const int, const int, const int);

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

        // Field内の値へのアクセスを wrap する
        tdArray& getScalar(const std::string varname) { return field->getScalar(varname); }

        void  setParent(Grid* g){ parent = g; }
        Grid* getParent(void){ return parent; }

        // Field 初期化
        void initializeField(void);

        //! 物体定義初期化
        void initializeObject(void);

        // 親子でのScalarやりとり用
        void copyScalarToChildren(std::string);
        void copyScalarToParent(std::string);

        // 子供管理メソッド
        void makeChild(const int, const int, const int, const int, const int, const int);
        void addChild(Grid*);
        void removeChild(const int);
        void checkGridValidness(void);

        std::vector<Grid*>& getChildren(void) {
            // 参照にしないと新しいポインタが生まれてしまう？
            return children;
        }

        int getChildrenLength(void) const {
            return children.size();
        }

        // child gridの総数を保存するためのメソッド
        int getSumOfChild(void) const {
            return sumTotalNumOfChildGrids;
        }

        void setSumOfChild(const int s) {
            sumTotalNumOfChildGrids = s;
        }

        void incrementSumOfChild(void);
        void decrementSumOfChild(void);

        //! Particleを格納したstd::vectorを、粒子種ごとに保持したstd::vector
        ParticleArray particles;

        // update fields
        void updateRho(void);

        void solvePoisson(void) {
            const int DEFAULT_ITERATION_LOOP = 500;
            field->solvePoisson(DEFAULT_ITERATION_LOOP, dx);
        }

        void updateEfield(void) {
            field->updateEfield(dx);
        }

        void updateEfieldFDTD(void) {
            field->updateEfieldFDTD(dx, dt);
        }

        void updateBfield(void) {
            field->updateBfield(dx, nx, ny, nz, dt);
        }

        // update particles
        void updateParticleVelocity(void);
        void updateParticlePosition(void);

        //! 粒子注入
        void injectParticles(void);

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

        // Extentを指定されたポインタに格納する (Silo MRG Tree出力用)
        void addExtent(int* logicalExtent[6], float* spatialExtent[6], float* rank[1]);

        // QuadMeshとVarをDBfileに突っ込む
        void putQuadMesh(DBfile*, std::string, const char* coordnames[3], int, DBoptlist*, DBoptlist*);

        // 出力用に密度を計算して返す
        float* getDensity(const int);

        // std out
        friend std::ostream& operator<<(std::ostream&, Grid*);

        //! inline functions
        void checkXBoundary(ParticleArray& pbuff, Particle& p, const double slx) {
            if(p.x < 0.0) {
                if (Environment::isPeriodic(AXIS::x, AXIS_SIDE::low)) {
                    // 計算空間の端でない場合は粒子を隣へ送る
                    // 計算空間の端にいるが、周期境界の場合も粒子を送る必要がある -> isPeriodic()でまとめて判定できる
                    pbuff[0].push_back(p);
                }
                p.isValid = 0;
            } else if (p.x > (slx - dx)) {
                if (Environment::isPeriodic(AXIS::x, AXIS_SIDE::up)) {
                    //! 計算空間の端でない場合、slx - dx から slx までの空間は下側の空間が担当するため、 slx を超えた場合のみ粒子を送信する
                    //! また、計算空間の端にいるが、周期境界の場合も同様の処理でよいため、isPeriodic()でまとめて判定できる
                    if (p.x > slx) {
                        pbuff[1].push_back(p);
                        p.isValid = 0;
                    }
                } else {
                    p.isValid = 0;
                }
            }
        }

        void checkYBoundary(ParticleArray& pbuff, Particle& p, const double sly) {
            if(p.y < 0.0) {
                if (Environment::isPeriodic(AXIS::y, AXIS_SIDE::low)) pbuff[2].push_back(p);
                p.isValid = 0;
            } else if (p.y > (sly - dx)) {
                if (Environment::isPeriodic(AXIS::y, AXIS_SIDE::up)) {
                    if (p.y > sly) {
                        pbuff[3].push_back(p);
                        p.isValid = 0;
                    }
                } else {
                    p.isValid = 0;
                }
            }
        }

        void checkZBoundary(ParticleArray& pbuff, Particle& p, const double slz) {
            if(p.z < 0.0) {
                if (Environment::isPeriodic(AXIS::z, AXIS_SIDE::low)) pbuff[4].push_back(p);
                p.isValid = 0;
            } else if (p.z > (slz - dx)) {
                if (Environment::isPeriodic(AXIS::z, AXIS_SIDE::up)) {
                    if (p.z > slz) {
                        pbuff[5].push_back(p);
                        p.isValid = 0;
                    }
                } else {
                    p.isValid = 0;
                }
            }
        }
};
#endif
