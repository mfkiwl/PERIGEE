#include "IonicModel_Test.hpp"

IonicModel_Test::IonicModel_Test()
  : IonicModel(0.1, 0.0, 140.0, 1.0)
{
  //SYS_T::commPrint("Test constructor. \n");
};

IonicModel_Test::~IonicModel_Test()
{
  //SYS_T::commPrint("Test destructor. \n");
};

void IonicModel_Test::print_info () const
{
  PetscPrintf(PETSC_COMM_WORLD,
	      "\t  Test model, no ionic current \n");
};

void IonicModel_Test::get_Iion(const double &r_old_in,
			      const double &V_in,
			      const double &I_stim,
			      double &f_r,
			      double &Iion) const
{

  Iion = 0.0 + I_stim/chi;
  f_r  = 0.0;
}


// EOF
