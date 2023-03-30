#ifndef _DLL_H_
#define _DLL_H_

#if BUILDING_DLL
# define DLLIMPORT __declspec (dllexport)
#else /* Not BUILDING_DLL */
# define DLLIMPORT __declspec (dllimport)
#endif /* Not BUILDING_DLL */

#define pi 3.1415927
#define COSZERO 1-(1e-12)
#define COS90D 1e-6
#define WTH 1e-4
#define CHANCE 0.1

struct SettingStruct{
       int np; //number of photons
       double dz; //delta z
       double dr; //delta r
       double da; //delta alpha
       int nz; //number of dz
       int nr; //number of dr
       int na; //number of da
       int nl; //number of layer
       double n_above; //refractive index of ambian above
       double n_below; //refractive index of ambian below
};
/*-------------------------------------------------*/
struct PhotonStruct //struct for a single photon
{
       double weight; //weight of photon
       int dead; //1 if dead or 0 if not
       int current_layer; //current layer
       double step_size; //step size
       double step_left; //step left
       double x; //x coordinate
       double y; //y coordinate
       double z; //z coordinate
       double ux; //directional cosine of x axis
       double uy; //directional cosine of y axis
       double uz; //directional cosine of z axis
       };




DLLIMPORT void MC(struct SettingStruct *setting,double *n,double *mua,double *mus,double *g,double *d,double *rsp,double *rd_ra,double *ab_rz,double *tt_ra);

DLLIMPORT double MC_Random();

DLLIMPORT void MC_SetElement(double *matrix,int indx,int indy,int xdim, double value);

DLLIMPORT double MC_GetElement(double *matrix,int indx,int indy,int xdim);

DLLIMPORT double *MC_SetCrit(double *n,struct SettingStruct setting);

DLLIMPORT double MC_Rspecular (double n0, double n1);

DLLIMPORT void MC_LaunchPhoton(struct PhotonStruct *photonPtr, double *rsp);

DLLIMPORT void MC_StepSize(struct PhotonStruct *photonPtr,double *muaPtr,double *musPtr);

DLLIMPORT int MC_HitBoundary(struct PhotonStruct *photonPtr,double *muaPtr,double *musPtr,double *dPtr);

DLLIMPORT void MC_MoveS(struct PhotonStruct *photonPtr);

DLLIMPORT void MC_Absorb(struct PhotonStruct *photonPtr,struct SettingStruct *setting,double *muaPtr,double *musPtr,double *ab_rz);

DLLIMPORT void MC_Scatter(struct PhotonStruct *photonPtr,double *gPtr);

DLLIMPORT double MC_RFresnel(double n1,double n2,double ca1,double *ca2_Ptr);

DLLIMPORT void MC_RecordR(double ref,struct PhotonStruct *photonPtr,struct SettingStruct *setting,double *rd_ra);

DLLIMPORT void MC_RecordT(double ref,struct PhotonStruct *photonPtr,struct SettingStruct *setting,double *tt_ra);

DLLIMPORT void MC_CrossUpOrNot(struct PhotonStruct *photonPtr,struct SettingStruct *setting,double *n,double *crit,double *rd_ra);

DLLIMPORT void MC_Roulette(struct PhotonStruct *photonPtr);

#endif /* _DLL_H_ */



















