#include "FIELD_3D.h"

FIELD_3D::FIELD_3D() :
  _xRes(0), _yRes(0), _zRes(0), _data(NULL)
{
}

FIELD_3D::FIELD_3D(int xRes, int yRes, int zRes) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _cellCount(xRes*yRes*zRes)
{
  _data = new particleVector[_xRes * _yRes * _zRes];
}

FIELD_3D::~FIELD_3D()
{
  if (_data) delete[] _data;
}

