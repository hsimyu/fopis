#ifndef __TDPIC_PARTICLE_H_INCLUDED__
#define __TDPIC_PARTICLE_H_INCLUDED__
#include <vector>
#include <string>

class Particle;

class Velocity {
    private:
        double vx, vy, vz;
    public:
        Velocity();
        ~Velocity();
        void set(double, double, double);
        double getVX(void) const;
        double getVY(void) const;
        double getVZ(void) const;
};

class Position {
    private:
        void updateDelta(void);

    public:
        int i, j, k;
        double x, y, z;
        double dx1, dy1, dz1, dx2, dy2, dz2;
        Position(const double, const double, const double);
        Position(const int, const int, const int);
        Position(const Particle&);
        ~Position();
        void setXYZ(const double, const double, const double);
        void setIJK(const int, const int, const int);
};

class Particle {
    public:
        // MPIで送るためにoffsetofを使う必要があり、
        // publicにする必要がある
        int typeId;
        int isValid;
        double x, y, z;
        double vx, vy, vz;

        Particle();
        Particle(Particle const&);
        ~Particle();

        void setPosition(Position const&);
        void setPosition(const double, const double, const double);

        void setVelocity(Velocity const&);
        void setVelocity(const double, const double, const double);

        double getEnergy(void) const;
        double getSquaredMagnitudeOfVelocity(void) const;

        void updatePosition(void);
        Particle& operator=(Particle const&);
        friend std::ostream& operator<<(std::ostream&, Particle const&);
};

class ParticleType {
    private:
        int id;
        std::string name;

        // plasma info
        std::string type;
        double charge;
        double mass;
        double density;
        double temperature;
        int size;

        // for computation
        int particle_per_cell;
        int totalNumber;
    public:
        ParticleType();

        void setId(int);
        void setCharge(double);
        void setMass(double);
        void setDensity(double);
        void setTemperature(double);
        void setSize(int);
        void setName(std::string);
        void setType(std::string);
        void setTotalNumber(int);
        void setPcell(int);

        int getId() const;
        std::string getType() const;
        std::string getName() const;
        double getCharge() const;
        double getMass() const;
        double getDensity() const;
        double getTemperature() const;
        double getTrueTemperature() const;
        int getSize() const;
        int getTotalNumber() const;
        int getPcell() const;

        int calcSize(void);
        int calcTotalNumber(void);
        double calcDeviation(void) const;
        std::string calcMemory(void) const;

        friend std::ostream& operator<<(std::ostream&, const ParticleType&);
        friend std::istream& operator<<(std::istream&, const ParticleType&);
};

typedef std::vector< std::vector<Particle> > ParticleArray;
#endif
