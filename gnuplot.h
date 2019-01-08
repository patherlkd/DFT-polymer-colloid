#ifndef GNU_H
#define GNU_H

class gnuplot
{ 

 public:
  gnuplot();
  ~gnuplot();
  void makefile(string);
  template <class T>
    void send2file(T x,T y);
  template <class T>
    void send3file(T x,T y,T z); 
  void plot(char);
  void splot(int);
  void xrange(db,db);
  void yrange(db,db);
  void ylabel(string);
  void xlabel(string);
  void erase();
 private:
  void print(string say){fprintf(gp,say.c_str());}
  FILE *gp;
  string file;
  ofstream out;
};

void gnuplot::xrange(db low,db high)
{ 
  fprintf(gp,"set xrange [%lf:%lf]\n",low,high);
} 


void gnuplot::yrange(db low,db high)
{ 
  fprintf(gp,"set yrange [%lf:%lf]\n",low,high);
} 


void gnuplot::xlabel(string say)
{ 
  fprintf(gp,"set xlabel '%s'\n",say.c_str());
} 


void gnuplot::ylabel(string say)
{ 
  fprintf(gp,"set ylabel '%s'\n",say.c_str());
} 


void gnuplot::erase() 
{
  fprintf(gp,"set key off\n"); 
  fprintf(gp,"set yrange [-1:1]\n");
  fprintf(gp,"plot 1/0\n"); 
} 
void gnuplot::splot(int pt)
{
  fprintf(gp, "splot '%s' u 1:2:3 w p pt %d lc rgb 'black'\n",file.c_str(),pt);

}
void gnuplot::plot(char mode)
{

  fprintf(gp,"plot '%s' u 1:2 w %c lc rgb 'black' \n",file.c_str(),mode);
  
  fflush(gp);
 
}

template <class T>
void gnuplot::send3file(T x,T y,T z)
{
  out << x <<'\t'<<y<<'\t'<<z<<endl; 
}
template <class T>
void gnuplot::send2file(T x, T y)
{
  out << x << '\t' << y <<endl;
}

void gnuplot::makefile(string name)
{
  file = name;
  out.open(file.c_str());
}

gnuplot::gnuplot()
{
FILE *p=popen("gnuplot -persist","w");
 gp=p;
}

gnuplot::~gnuplot() 
{ 
  fclose(gp);
}


#endif
