#include "jetBackgroundCut.h"
#include <ffaobjects/EventHeaderv1.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfoContainerSimv1.h>
#include <calobase/TowerInfoContainerSimv2.h>
#include <globalvertex/GlobalVertexMapv1.h>
#include <globalvertex/GlobalVertex.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtContainerV1.h>
#include <mbd/MbdPmtSimContainerV1.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <globalvertex/MbdVertexMapv1.h>
#include <globalvertex/MbdVertex.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <ffaobjects/EventHeaderv1.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <HepMC/GenEvent.h>
#include <mbd/MbdPmtHit.h>
#include <jetbackground/TowerBackgroundv1.h>
#include <cmath>
#include <mbd/MbdOut.h>
#include <TLorentzVector.h>
#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>
#include <g4centrality/PHG4CentralityReco.h>
#include <jetbase/JetContainerv1.h>
#include <jetbase/Jet.h>
#include <calobase/RawTowerv1.h>
#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/SubtractTowersCS.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <iostream>

#include <centrality/CentralityInfo.h>
#include <calotrigger/MinimumBiasInfo.h>
#include <calotrigger/MinimumBiasInfov1.h>
#include <ffarawobjects/Gl1Packetv2.h>
#include <jetbase/JetMapv1.h>
#include <jetbase/JetMap.h>
using namespace std;
//____________________________________________________________________________..
jetBackgroundCut::jetBackgroundCut(const std::string &name, const int debug, int datorsim, int dotow):
  SubsysReco("test")//).c_str())
{
  _dotow = dotow;
  _evtnum = 0;
  _evtct = 0;
  _foutname = name;
  _debug = debug;
  mbevt = 0;
  _datorsim = datorsim;
}

//____________________________________________________________________________..
jetBackgroundCut::~jetBackgroundCut()
{

}

//____________________________________________________________________________..
int jetBackgroundCut::Init(PHCompositeNode *topNode)
{

  _rc = recoConsts::instance();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int jetBackgroundCut::process_event(PHCompositeNode *topNode)
{

  _cutBits = 0;
  TowerInfoContainer *towersEM = findNode::getClass<TowerInfoContainerSimv1>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  TowerInfoContainer *towersIH = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainer *towersOH = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_HCALOUT");
  JetContainer *jets = findNode::getClass<JetContainerv1>(topNode, "AntiKt_Tower_HIRecoSeedsRaw_r04");
  MbdVertexMap* mbdvtxmap = findNode::getClass<MbdVertexMapv1>(topNode, "MbdVertexMap");
  GlobalVertexMap* gvtxmap = findNode::getClass<GlobalVertexMapv1>(topNode, "GlobalVertexMap");

  RawTowerGeomContainer *geom[3];
  geom[0] = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_CEMC");
  geom[1] = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_HCALIN");
  geom[2] = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_HCALOUT");


  vtx[0] = 0;
  vtx[1] = 0;
  vtx[2] = NAN;
  if(_datorsim || !gvtxmap)
    {
      for(auto iter = mbdvtxmap->begin(); iter != mbdvtxmap->end(); ++iter)
	{
	  MbdVertex* mbdvtx = iter->second;
	  vtx[2] = mbdvtx->get_z();
	  break;
	}
    }
  else if(gvtxmap)
    {
      auto iter = gvtxmap->begin();
      while(iter != gvtxmap->end())
	{
	  GlobalVertex* gvtx = iter->second;
	  vtx[2] = gvtx->get_z();
	  vtx[0] = gvtx->get_x();
	  vtx[1] = gvtx->get_y();
	  iter++;
	  break;
	}
    }
  if(std::isnan(vtx[2]) || abs(vtx[2]) > 150)
    {
      if(ismb) mbevt--;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  if(_debug > 1) cout << "Getting jets: " << endl;
  if(jets)
    {
      int tocheck = jets->size();
      if(_debug > 2) cout << "Found " << tocheck << " jets to check..." << endl;
      for(int i=0; i<tocheck; ++i)
	{
	  Jet *jet = jets->get_jet(i);
	  if(jet)
	    {
	      jet_et[njet] = jet->get_eta();
	      jet_e[njet] = jet->get_e()/cosh(jet_et[njet]);
	    }
	  else
	    {
	      continue;
	    }
	  _frcem[njet] = 0;
	  _frcoh[njet] = 0;
	  if(jet_e[njet] < 8) continue;
	  if(_debug > 2) cout << "found a good jet!" << endl;
	  float maxeovertot = 0;
	  float hcale = 0;
	  float ihcale = 0;
	  float ohcale = 0;
	  float ecale = 0;
	  int ncomp = 0;
	  
	  for(auto comp: jet->get_comp_vec())
	    {
	      unsigned int channel = comp.second;
	      TowerInfo* tower;
	      if(comp.first == 13 || comp.first == 28 || comp.first == 25)
		{
		  tower = towersEM->get_tower_at_channel(channel);
		  int key = towersEM->encode_key(channel);
		  const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, towersEM->getTowerEtaBin(key), towersEM->getTowerPhiBin(key));
		  RawTowerGeom *tower_geom = geom[1]->get_tower_geometry(geomkey);
		  float radius = 93.5;
		  float ihEta = tower_geom->get_eta();
		  float emZ = radius/(tan(2*atan(exp(-ihEta))));
		  float newz = emZ - vtx[2];
		  float newTheta = atan2(radius,newz);
		  float towerEta = -log(tan(0.5*newTheta));
		  _frcem[njet] += tower->get_energy()/cosh(towerEta);
		}
	      if(comp.first == 7 || comp.first == 27)
		{
		  tower = towersOH->get_tower_at_channel(channel);
		  int key = towersOH->encode_key(channel);
		  const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, towersOH->getTowerEtaBin(key), towersOH->getTowerPhiBin(key));
		  RawTowerGeom *tower_geom = geom[2]->get_tower_geometry(geomkey);
		  float radius = tower_geom->get_center_radius();
		  float newz = tower_geom->get_center_z() - vtx[2];
		  float newTheta = atan2(radius,newz);
		  float towerEta = -log(tan(0.5*newTheta));
		  _frcoh[njet] += tower->get_energy()/cosh(towerEta);
		}
	    }
	  _frcem[njet] /= jet_e[njet];
	  _frcoh[njet] /= jet_e[njet];

  else
    {
      if(_debug > 1) cout << "No jet node!" << endl;
    }
  if(!njet && !_dotow) return Fun4AllReturnCodes::EVENT_OK;
  
  return Fun4AllReturnCodes::EVENT_OK;
    
}
//____________________________________________________________________________..
int jetBackgroundCut::ResetEvent(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "jetBackgroundCut::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int jetBackgroundCut::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "jetBackgroundCut::End(PHCompositeNode *topNode) This is the End..." << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int jetBackgroundCut::Reset(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "jetBackgroundCut::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void jetBackgroundCut::Print(const std::string &what) const
{
  std::cout << "jetBackgroundCut::Print(const std::string &what) const Printing info for " << what << std::endl;
}
