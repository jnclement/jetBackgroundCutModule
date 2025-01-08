#ifndef JETBACKGROUNDCUT_H
#define JETBACKGROUNDCUT_H

#include <fun4all/SubsysReco.h>
#include <string>
#include <cmath>
#include <phool/recoConsts.h>

class PHCompositeNode;
class CentralityInfo;
class jetBackgroundCut : public SubsysReco
{
 public:

  jetBackgroundCut(const std::string &name = "jetBackgroundCutModule", const bool debug = 0, const bool doAbort = 0);

  virtual ~jetBackgroundCut();

  bool failsLoEmFracETCut(float emFrac, float ET, bool dPhiCut, bool isDijet)
  {
    return (frcem < 0.1 && jet_ET > (50*frcem+20)) && (dPhiCut || !isDijet);
  }

  bool failsHiEmFracETCut(float emFrac, float ET, bool dPhiCut, bool isDijet);
  {
    return (frcem > 0.9 && jet_ET > (-50*frcem+70)) && (dPhiCut || !isDijet);
  }

  bool failsIhFracCut(float emFrac, float ohFrac);
  {
    return frcem + frcoh < 0.65;
  }

  bool failsdPhiCut(float dPhi, bool isDijet)
  {
    return dPhi < 3*M_PI/4 && isDijet;
  }

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;


 private:
  recoConsts *_rc;
  bool _doAbort;
  string _name;
  bool _debug;
};

#endif // R24TREEMAKER
