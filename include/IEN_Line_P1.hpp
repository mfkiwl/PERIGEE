#ifndef IEN_LINE_P1_HPP
#define IEN_LINE_P1_HPP
// ==================================================================
// IEN_Line_P1.hpp
//
// This objects gives the IEN array for a 3D linear line
// mesh.
//
// This mesh is assumed to be linear line elements. Hence, the 
// number of local nodes is 2. The length of the IEN array is 2 nElem.
//
// Author: Ju Liu
// Modified: Oguz Ziya Tikenogullari - from the original: IEN_Tetra_P1
// Date: Dec.18 2016.
// ==================================================================
#include <vector>
#include "IIEN.hpp"

class IEN_Line_P1 : public IIEN
{
  public:
    // Constructor: assume the ien_array is read from the .vtu file
    IEN_Line_P1( const int &in_nelem, const std::vector<int> &in_ien );
    
    ~IEN_Line_P1();

    virtual int get_IEN( const int &ee, const int &l_node ) const;

    virtual void print_IEN() const;

  private:
    int * IEN;

    int nElem;
};

#endif
