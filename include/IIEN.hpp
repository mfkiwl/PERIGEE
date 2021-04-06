#ifndef IIEN_HPP
#define IIEN_HPP
// ==================================================================
// IIEN.hpp
// The interface for IEN array classes.
//
// Date: Sept. 24th 2013
// ==================================================================
#include <cstdlib>
#include <iostream>

class IIEN
{
  public:
    IIEN(){};
    virtual ~IIEN(){};

    // get the IEN arrray for element e at local node l_node
    virtual int get_IEN( const int &ee, const int &l_node ) const = 0;

    // print IEN array
    virtual void print_IEN() const = 0;

    //get total number of functions, useful for meshes of non-uniform  elemtype
    virtual int get_nFunc_tot() const
    {std::cerr<<"Error: get_nFunc_tot is not implemented. \n"; exit(EXIT_FAILURE);}      

    //get total number of elements, useful for meshes of non-uniform  elemtype
    virtual int get_nElem_tot() const
    {std::cerr<<"Error: get_nElem_tot is not implemented. \n"; exit(EXIT_FAILURE);}      
		
    // print info
    virtual void print_info() const
    {std::cerr<<"Error: print_info is not implemented. \n"; exit(EXIT_FAILURE);}

    // get the IEN pointer
    virtual int * get_IENptr() const
    {std::cerr<<"Error: get_IENptr is not implemented. \n"; exit(EXIT_FAILURE); return NULL;}
};

#endif
