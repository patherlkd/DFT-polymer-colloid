#ifndef USEFUL_H
#define USEFUL_H



typedef double db;
typedef long double dbl;

//FUNCTIONS

// CONSTANTS
const db pi = 3.14159265359; 
const db euler = 2.71828; 
//useful CLASSES


//{} ARRAY
class dbarray
{
public:
  dbarray();
  ~dbarray();
  void setelem(int,db);
  void printall();
  void resize(int);
  void setnumelems(int);
  void zero();
  void addelem(int);
  void divideelems(db);
  db getelem(int);
private:

  db *elems;
  size_t numelems;

};

#endif
