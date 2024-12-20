#ifndef JETBACKGROUNDCUT_H
#define JETBACKGROUNDCUT_H

#include <fun4all/SubsysReco.h>
#include <string>
#include <phool/recoConsts.h>

class PHCompositeNode;
class CentralityInfo;
class jetBackgroundCut : public SubsysReco
{
 public:

  jetBackgroundCut(const std::string &name = "jetBackgroundCut", const int debug = 0, const int doAbort = 0);

  virtual ~jetBackgroundCut();

  int failsEmFracHiETCut(float emFrac, float ET);

  int failsEmFracLoETCut(float emFrac, float ET);
 
  int failsIhFracCut(float emFrac, float ohFrac);

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;


 private:
  unsigned char _cutBits;
  _recoConsts *rc;
};

#endif // R24TREEMAKER
