#include "FIELD_3D.h"

FIELD_3D::FIELD_3D() :
  _xRes(0), _yRes(0), _zRes(0), _data(NULL)
{
}

FIELD_3D::FIELD_3D(int xRes, int yRes, int zRes) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _cellCount(xRes*yRes*zRes)
{
  _data = new particleVector[_xRes * _yRes * _zRes];
  for (int i = 0; i < _cellCount; i++)
    _data[i].reserve(1000);
  // clear();
  cout << "grid is " << xRes << "x" << yRes << "x" << zRes << endl;
}

FIELD_3D::~FIELD_3D()
{
  //cout << "field deleted" << endl;
  if (_data) delete[] _data;
  // delete the actual vectors?
}

