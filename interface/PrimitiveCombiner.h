#ifndef L1Trigger_L1IntegratedMuonTrigger_PrimitiveCombiner_h_
#define L1Trigger_L1IntegratedMuonTrigger_PrimitiveCombiner_h_
// 
// Class: L1ITMu:: PrimitiveCombiner
//
// Info: This class combine information from DT and/or RPC primitives
//       in order to calculate better phi/phiBending
//
// Author: 
//
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include <cmath>

namespace L1ITMu {

  class TriggerPrimitive;

  class PrimitiveCombiner {

  public :
    /// a struct useful for resulution info sharing
    struct resolutions {
      double xDt;
      double xRpc;
      double phibDtCorr;
      double phibDtUnCorr;
      resolutions( const double & xDtResol, const double & xRpcResol,
		   const double & phibDtCorrResol, const double & phibDtUnCorrResol )
	:  xDt( xDtResol ), xRpc( xRpcResol ),
	   phibDtCorr( phibDtCorrResol ), phibDtUnCorr( phibDtUnCorrResol ) {};
    };

  public :
    explicit PrimitiveCombiner( const DTChamberId & dtId, 
				const resolutions & res );

    /// feed the combiner with the available primitives
    void addDt( const L1ITMu::TriggerPrimitive & prim );
    void addDtHI( const L1ITMu::TriggerPrimitive & prim );
    void addDtHO( const L1ITMu::TriggerPrimitive & prim );
    void addRpcIn( const L1ITMu::TriggerPrimitive & prim );
    void addRpcOut( const L1ITMu::TriggerPrimitive & prim );

    /// do combine the primitives
    void combine();

    /// output result variables
    int bx() const { return _bx;};
    int radialAngle() const { return _radialAngle;};
    int bendingAngle() const { return _bendingAngle;};
    bool isValid() const { return _isValid;};


    /// FIXME : Calculates new phiBending, check how to use 
    inline float phiBCombined( const float & xDt, const float & zDt,
			       const float & xRpc, const float & zRpc )
      {
	return (xRpc - xDt) / (zRpc - zDt);
      };
    /// FIXME END

    /// FIXME : Calculates new phiBending resolution
    inline float phiBCombinedResol( const float & resol_xDt,
				    const float & resol_xRpc )
      {
	return sqrt( resol_xRpc*resol_xRpc + resol_xDt*resol_xDt );
      };
    /// FIXME END

  private :

    /// a struct for internal usage: store results
    struct results {
      double radialAngle;
      double bendingAngle;
      double bendingResol;
      double weight;
      results() : radialAngle(0), bendingAngle(0), bendingResol(0), weight(0) {};
    };


    /// Calculates new phiBending, check how to use weights
    results combineDt( const L1ITMu::TriggerPrimitive * dt,
		       const L1ITMu::TriggerPrimitive * rpc );

    /// Calculates new phiBending, check how to use weights
    results combineDtRpc( const L1ITMu::TriggerPrimitive * dt,
			  const L1ITMu::TriggerPrimitive * rpc );

  private :
    DTChamberId _dtId;
    resolutions _resol;

    int _bx;
    int _radialAngle;
    int _bendingAngle;
    int _bendingResol;
    bool _isValid;

    const L1ITMu::TriggerPrimitive * _dtHI;
    const L1ITMu::TriggerPrimitive * _dtHO;
    const L1ITMu::TriggerPrimitive * _rpcIn;
    const L1ITMu::TriggerPrimitive * _rpcOut;

  };
}

#endif
