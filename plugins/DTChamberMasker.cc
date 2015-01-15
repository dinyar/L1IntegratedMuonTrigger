// -*- C++ -*-
//
// Package:    L1Trigger/DTChamberMasker
// Class:      DTChamberMasker
// 
/**\class DTChamberMasker DTChamberMasker.cc L1Trigger/DTChamberMasker/plugins/DTChamberMasker.cc

 Description: Class to mask DT digis or trigger segments on a chamber by chamber basis

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Carlo Battilana
//         Created:  Sun, 11 Jan 2015 15:12:51 GMT
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <regex>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"

#include "Geometry/DTGeometry/interface/DTChamber.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"

#include "DataFormats/DTDigi/interface/DTLocalTriggerCollection.h"
#include "DataFormats/DTDigi/interface/DTDigiCollection.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"

//
// class declaration
//

class DTChamberMasker : public edm::EDProducer {
public:
  explicit DTChamberMasker(const edm::ParameterSet&);
  ~DTChamberMasker();

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
      
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;

  void createMaskedChamberCollection(edm::ESHandle<DTGeometry> &); 

  // ----------member data ---------------------------

  edm::InputTag m_digiTag;
  edm::InputTag m_trPhTag;
  edm::InputTag m_trThTag;

  bool m_doTriggerFromDDU;
  
  std::vector<std::regex> m_maskedChRegEx;
  std::vector<DTChamberId> m_enabledCh;  
  
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
DTChamberMasker::DTChamberMasker(const edm::ParameterSet& iConfig) :
  m_digiTag( iConfig.getParameter<edm::InputTag>("digiTag") ),
  m_trPhTag( iConfig.getParameter<edm::InputTag>("triggerPrimPhTag") ),
  m_trThTag( iConfig.getParameter<edm::InputTag>("triggerPrimThTag") ),
  m_doTriggerFromDDU( iConfig.getParameter<bool>("doTriggerFromDDU") ) 
{

  produces<DTDigiCollection>();
  if (m_doTriggerFromDDU)
    produces<DTLocalTriggerCollection>();
  
  produces<L1MuDTChambPhContainer>();
  produces<L1MuDTChambThContainer>();

  for ( auto regStr : iConfig.getParameter<std::vector<std::string>>("maskedChRegEx") )
    m_maskedChRegEx.push_back(std::regex(regStr));

}


DTChamberMasker::~DTChamberMasker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DTChamberMasker::produce(edm::Event& event, const edm::EventSetup& config)
{

  // CB Fix a few things:
  // add code to mask DDU trigger primitives as well
  // add possibility to put hit efficiency at a different than 100 or 0
  // deal better with SEC 13 and 14 for trigger
  // change enabled object maps from ID to rawId ?
  // add FillDescriptors
  // write loops better
  
  std::auto_ptr<DTDigiCollection> filteredDigis(new DTDigiCollection());
  std::auto_ptr<DTLocalTriggerCollection> filteredTrigsDDU(new DTLocalTriggerCollection());

  std::auto_ptr<L1MuDTChambPhContainer> filteredPhTrigPrimsCont(new L1MuDTChambPhContainer());
  std::auto_ptr<L1MuDTChambThContainer> filteredThTrigPrimsCont(new L1MuDTChambThContainer());

  std::vector<L1MuDTChambPhDigi> filteredPhTrigPrims;
  std::vector<L1MuDTChambThDigi> filteredThTrigPrims;
  

  if (m_digiTag.label() != "")
    {
  
      edm::Handle<DTDigiCollection> dtDigis;
      event.getByLabel(m_digiTag, dtDigis);
  
      DTDigiCollection::DigiRangeIterator dtLayerIdIt  = dtDigis->begin();
      DTDigiCollection::DigiRangeIterator dtLayerIdEnd = dtDigis->end();
      
      for (; dtLayerIdIt != dtLayerIdEnd; ++dtLayerIdIt)
	{
	  
	  DTChamberId chId = ((*dtLayerIdIt).first).chamberId();
	  if (std::find(m_enabledCh.begin(),m_enabledCh.end(),chId) != m_enabledCh.end())
	    filteredDigis->put((*dtLayerIdIt).second,(*dtLayerIdIt).first);
	  
	}
      // CB write the code to apply masks to the DDU collection here
      //      if (m_doTriggerFromDDU)
      //	{
      //        }

    }

  if (m_trPhTag.label() != "")
    {
  
      edm::Handle<L1MuDTChambPhContainer> phTrigPrimDigis;
      event.getByLabel(m_trPhTag, phTrigPrimDigis);
      
      for (auto const trPhDigi : (*phTrigPrimDigis->getContainer()) )
	{
	  
	  DTChamberId chId = DTChamberId(trPhDigi.whNum(), trPhDigi.stNum(), (trPhDigi.scNum()+1));
	  if (std::find(m_enabledCh.begin(),m_enabledCh.end(),chId) != m_enabledCh.end())
	    filteredPhTrigPrims.push_back(trPhDigi);	  
	  
	}

    }
      
  if (m_trThTag.label() != "")
    {
  
      edm::Handle<L1MuDTChambThContainer> thTrigPrimDigis;
      event.getByLabel(m_trThTag, thTrigPrimDigis);

      for (auto const trThDigi : (*thTrigPrimDigis->getContainer()) )
	{
	  
	  DTChamberId chId = DTChamberId(trThDigi.whNum(), trThDigi.stNum(), (trThDigi.scNum()+1));
	  if (std::find(m_enabledCh.begin(),m_enabledCh.end(),chId) != m_enabledCh.end())
	    filteredThTrigPrims.push_back(trThDigi);	  
	  
	}

    }

  filteredPhTrigPrimsCont->setContainer(filteredPhTrigPrims);
  filteredThTrigPrimsCont->setContainer(filteredThTrigPrims);
  
  event.put(filteredDigis);
  if (m_doTriggerFromDDU)
    event.put(filteredTrigsDDU);
  event.put(filteredPhTrigPrimsCont);
  event.put(filteredThTrigPrimsCont);

}

// ------------ method called once each job just before starting event loop  ------------
void 
DTChamberMasker::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
DTChamberMasker::endJob()
{

}

// ------------ method called when starting to processes a run  ------------
void
DTChamberMasker::beginRun(edm::Run const& run, edm::EventSetup const& config)
{

  edm::ESHandle<DTGeometry> dtGeom;
  config.get<MuonGeometryRecord>().get(dtGeom);

  if (m_enabledCh.size() == 0)
    {

      createMaskedChamberCollection(dtGeom);

      std::cout << "[DTChamberMasker::beginRun]: List of enabled chambers: " << std::endl;

      for ( const DTChamberId chId : m_enabledCh)
	std::cout << "[DTChamberMasker::beginRun]:" << chId << std::endl;

    }

}

  

void
DTChamberMasker::createMaskedChamberCollection(edm::ESHandle<DTGeometry> & dtGeom)
{

  const std::vector<const DTChamber*>& chambers = dtGeom->chambers();

  for ( const DTChamber *ch : chambers)
   {

     DTChamberId chId = ch->id();

     std::string chTag = "WH" + std::to_string(chId.wheel())
                       + "_ST" + std::to_string(chId.station())
                       + "_SEC" + std::to_string(chId.sector());

     bool isMaskedCh = false;

     for ( const std::regex chRegex : m_maskedChRegEx)
       if ( std::regex_match(chTag,chRegex) )
	 isMaskedCh = true;

     if (!isMaskedCh)
       	 m_enabledCh.push_back(chId);

   }
  
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DTChamberMasker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DTChamberMasker);
