#include "IEN_Line_P1.hpp"

IEN_Line_P1::IEN_Line_P1( const int &in_nelem,
    const std::vector<int> &in_ien )
{
  nElem = in_nelem;

  // Linear line element has two local nodes
  IEN = new int [nElem * 2];

  for(unsigned int ii=0; ii<in_ien.size(); ++ii) IEN[ii] = in_ien[ii];
}


IEN_Line_P1::~IEN_Line_P1()
{
  delete [] IEN; IEN = NULL;
}


int IEN_Line_P1::get_IEN( const int &ee, const int &l_node ) const
{
  return IEN[ee*2+l_node];
}


void IEN_Line_P1::print_IEN() const
{
  std::cout<<std::endl;
  std::cout<<"====== IEN ====== \n";
  for(int ii=0; ii<nElem; ++ii)
  {
    for(int jj=0; jj<2; ++jj)
      std::cout<<get_IEN(ii, jj)<<'\t';
    std::cout<<std::endl;
  }
  std::cout<<"================= \n";
}

// EOF
