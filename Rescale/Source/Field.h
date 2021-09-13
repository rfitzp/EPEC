// Field.h

// ################
// pFile data class
// ################

#ifndef FIELD
#define FIELD

// ############
// Class header
// ############
class Field
{
 private:

  // ----------
  // Class data
  // ----------
  int     N;     // Number of data points
  double* X;     // Array of X values
  double* Y;     // Array of Y values
  double* dYdX;  // Array of dY/dX values
  
  // ----------------------
  // Public class functions
  // ----------------------
 public:
    
   Field ();      
   Field (int n); 
  ~Field ();      

  void resize (int n);
  
  int    GetN ();
  void   PushData (int i, double  x, double  y, double  dydx);
  void   PullData (int i, double& x, double& y, double& dydx);
  double GetX     (int i);
  double GetY     (int i);
  double GetdYdX  (int i);
  void   Rescale  (double A1);
  
  // -----------------------
  // Private class functions
  // -----------------------
 private:
};

#endif //FIELD
