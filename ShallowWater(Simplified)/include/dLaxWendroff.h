#ifndef DLAXWENDROFF_H
#define DLAXWENDROFF_H
//-----------------------------//
class dLaxWendroff {
public:
    explicit dLaxWendroff(unsigned int EquationsNumP);
    ~dLaxWendroff() = default;
private:
    unsigned int EquationsNum = 1;

    double TimeStep = 0.0;

//    virtual
};
//-----------------------------//
#endif
