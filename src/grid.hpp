#ifndef __TDPIC_GRID_H_INCLUDED__
#define __TDPIC_GRID_H_INCLUDED__
#include <vector>
#include <map>
#include <silo.h>
#include "particle.hpp"
#include "environment.hpp"
#include "utils.hpp"

class Field;
class Particle;

//! @class Grid
class Grid {
    private:
        Grid* parent;
        std::vector<Grid*> children;

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

        Field* field;

        // Class Unique ID
        static unsigned int nextID;

    public:
        Grid(void);
        Grid(Grid*, const int, const int, const int, const int, const int, const int);

        ~Grid();

        unsigned int getID(void) const;
        unsigned int getNextID(void);
        // IDへのsetterは提供しない、代わりにresetterを用意
        static void resetNextID(void);

        void setFromIX(int);
        void setFromIY(int);
        void setFromIZ(int);
        int  getFromIX() const;
        int  getFromIY() const;
        int  getFromIZ() const;
        void setToIX(int);
        void setToIY(int);
        void setToIZ(int);
        int  getToIX() const;
        int  getToIY() const;
        int  getToIZ() const;

        void setBaseX(double);
        void setBaseY(double);
        void setBaseZ(double);

        double getBaseX() const;
        double getBaseY() const;
        double getBaseZ() const;

        void setNX(int);
        void setNY(int);
        void setNZ(int);

        int  getNX() const;
        int  getNY() const;
        int  getNZ() const;

        void setLevel(int);
        int getLevel() const;

        void setDX(double);
        double getDX() const;

        void setParent(Grid*);
        Grid* getParent();

        void setField(Field*);
        Field* getField();

        // Field 初期化
        void initializeField(void);

        // 親子でのScalarやりとり用
        void copyScalarToChildren(std::string);
        void copyScalarToParent(std::string);

        // 子供管理メソッド
        void makeChild(const int, const int, const int, const int, const int, const int);
        void addChild(Grid*);
        std::vector<Grid*>& getChildren();
        int getChildrenLength(void) const;
        void removeChild(const int);
        void checkGridValidness(void);

        // child gridの総数を保存するためのメソッド
        int getSumOfChild(void) const;
        void setSumOfChild(const int);
        void incrementSumOfChild(void);
        void decrementSumOfChild(void);

        //! Particleを格納したstd::vectorを、粒子種ごとに保持したstd::vector
        ParticleArray particles;

        // update fields
        void updateRho(void);
        void solvePoisson(void);
        void updateEfield(void);
        void updateBfield(void);

        // update particles
        void updateParticleVelocity(void);
        void updateParticlePosition(void);

        //! エネルギーを計算する
        double getParticleEnergy(void);
        double getFieldEnergy(void);

        // create mesh nodes array
        float** getMeshNodes(int);
        int* getChildOfPatches(void);
        int* getNumOfPatches(void);
        int* getChildIdMap(void);
        int getMaxLevel(void);
        int getMaxChildLevel(void);

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

        // std out
        friend std::ostream& operator<<(std::ostream&, Grid*);

        //! inline functions
        void checkXBoundary(std::vector< std::vector<Particle> >& pbuff, Particle& p, const double slx) {
            if(p.x < 0.0) {
                if(!Environment::onLowXedge) pbuff[0].push_back(p);
                p.isValid = 0;
            } else if (p.x > slx) {
                if(!Environment::onHighXedge) pbuff[1].push_back(p);
                p.isValid = 0;
            }
        }

        void checkYBoundary(std::vector< std::vector<Particle> >& pbuff, Particle& p, const double sly) {
            if(p.y < 0.0) {
                if(!Environment::onLowYedge) pbuff[2].push_back(p);
                p.isValid = 0;
            } else if (p.y > sly) {
                if(!Environment::onHighYedge) pbuff[3].push_back(p);
                p.isValid = 0;
            }
        }

        void checkZBoundary(std::vector< std::vector<Particle> >& pbuff, Particle& p, const double slz) {
            if(p.z < 0.0) {
                if(!Environment::onLowZedge) pbuff[4].push_back(p);
                p.isValid = 0;
            } else if (p.z > slz) {
                if(!Environment::onHighZedge) pbuff[5].push_back(p);
                p.isValid = 0;
            }
        }
};
#endif
