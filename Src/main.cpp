# include "gflame.hpp"

#include <AMReX.H>

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    Gflame::init();
    amrex::Finalize();
}
