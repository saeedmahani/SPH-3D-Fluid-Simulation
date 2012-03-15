#include "FIELD_2D.h"

///////////////////////////////////////////////////////////////////////
// Constructor / Destructor
///////////////////////////////////////////////////////////////////////
FIELD_2D::FIELD_2D() : 
  _xRes(0), _yRes(0), _data(NULL) 
{ 
} 

FIELD_2D::FIELD_2D(int xRes, int yRes) : 
  _xRes(xRes), _yRes(yRes)
{
  _data = new double[_xRes * _yRes];
  clear();
}

FIELD_2D::~FIELD_2D() 
{ 
  if (_data) delete[] _data; 
}

void FIELD_2D::print() {
  for (int i = 0; i < _xRes * _yRes; i++) {
    printf("%.10f,", _data[i]);
  }
  printf("\n");
}

///////////////////////////////////////////////////////////////////////
// set the size of the array
///////////////////////////////////////////////////////////////////////
void FIELD_2D::resize(int xRes, int yRes) {
  if (_data) delete[] _data;
  _xRes = xRes;
  _yRes = yRes;
  _data = new double[_xRes * _yRes];
  clear();
}

///////////////////////////////////////////////////////////////////////
// set to zero
///////////////////////////////////////////////////////////////////////
void FIELD_2D::clear() {
  for (int x = 0; x < _xRes * _yRes; x++)
    _data[x] = 0.0; 
}

///////////////////////////////////////////////////////////////////////
// scale the field by a constant
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator*=(const double& scalar)
{
  for (int x = 0; x < _xRes * _yRes; x++)
    _data[x] *= scalar;

  return *this;
}

///////////////////////////////////////////////////////////////////////
// set one field to another
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator=(const FIELD_2D& field)
{
  assert(field.xRes() == _xRes);
  assert(field.yRes() == _yRes);

  for (int x = 0; x < _xRes * _yRes; x++)
    _data[x] = field._data[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
// squared sum of entries
///////////////////////////////////////////////////////////////////////
double FIELD_2D::squaredSum() const
{
  double final = 0.0;
  for (int x = 0; x < _xRes * _yRes; x++)
    final += _data[x] * _data[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
// dot product with another field
///////////////////////////////////////////////////////////////////////
double FIELD_2D::dotProduct(FIELD_2D& field) const
{
  assert(field.xRes() == _xRes);
  assert(field.yRes() == _yRes);

  double final = 0.0;
  for (int x = 0; x < _xRes * _yRes; x++)
    final += _data[x] * field._data[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
// subtract another field from this field
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator-=(const FIELD_2D& field)
{
  assert(field.xRes() == _xRes);
  assert(field.yRes() == _yRes);

  for (int x = 0; x < _xRes * _yRes; x++)
    _data[x] -= field._data[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
// add another field to this field
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator+=(const FIELD_2D& field)
{
  assert(field.xRes() == _xRes);
  assert(field.yRes() == _yRes);

  for (int x = 0; x < _xRes * _yRes; x++)
    _data[x] += field._data[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
// axpy operation with another field
// this = this + scalar * field
///////////////////////////////////////////////////////////////////////
void FIELD_2D::axpy(const double& scalar, const FIELD_2D& field)
{
  assert(field.xRes() == _xRes);
  assert(field.yRes() == _yRes);
  
  for (int x = 0; x < _xRes * _yRes; x++)
    _data[x] += scalar * field._data[x];
}
