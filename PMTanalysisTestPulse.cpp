#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TDirectory.h>
#include <TText.h>
#include <TLatex.h>
#include <TLine.h>
#include <TBox.h>
#include <TLegend.h>
#include <Riostream.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

//--------------------------------------------------------------------------------
// Structures
//--------------------------------------------------------------------------------
struct RunInfo {
  int runNumber;
  int width;
  std::string offsetStr;
  bool isOriginal;
  double numericOffset;
};

struct TowerDefault {
  int sector;
  int ib;
  int iphi;
  int ieta;
  double defaultOffset;
};

// Forward-declare the XYerr struct if needed:
struct XYerr {
    double x;
    double y;
    double e;
};

struct BoardFitData {
    bool isOld;
    int sector;
    int ib;

    // data[w] => vector<XYerr>
    std::map<int, std::vector<XYerr>> data;

    int bestW = -1;

    double amplitude    = 0.0;
    double amplitudeErr = 0.0;
    double kParam       = 0.0;
    double kParamErr    = 0.0;
};

// A helper struct if you want them declared here rather than inside doAnalysis:
struct Kval {
    int sector;
    int ib;
    double kOld;
    double kNew;
};
struct AmpVal {
    int sector;
    int ib;
    double ampOld; // A from old bestW fit
    double ampNew; // A from new bestW fit
};


//--------------------------------------------------------------------------------
// Mapping logic
//--------------------------------------------------------------------------------
int custom_sector_mapping(unsigned int eta, unsigned int phi)
{
  if(phi >= 256) return -1;
  int baseSector = phi / 8;
  int sector = -1;
  if(eta >= 48 && eta < 96) {
    sector = baseSector;
  } else if(eta < 48) {
    sector = 32 + baseSector;
  }
  if(sector < 0 || sector > 63) return -1;
  return sector;
}

int custom_ib_board(int eta, int phi)
{
  if(phi >= 256) return -1;
  int ib = -1;
  // top half
  if     (eta >= 48 && eta < 56) ib = 0;
  else if(eta >= 56 && eta < 64) ib = 1;
  else if(eta >= 64 && eta < 72) ib = 2;
  else if(eta >= 72 && eta < 80) ib = 3;
  else if(eta >= 80 && eta < 88) ib = 4;
  else if(eta >= 88 && eta < 96) ib = 5;
  // bottom half
  else if(eta >= 40 && eta < 48) ib = 0;
  else if(eta >= 32 && eta < 40) ib = 1;
  else if(eta >= 24 && eta < 32) ib = 2;
  else if(eta >= 16 && eta < 24) ib = 3;
  else if(eta >= 8  && eta < 16) ib = 4;
  else if(eta >= 0  && eta < 8 ) ib = 5;
  return ib;
}

//--------------------------------------------------------------------------------
// Misc utilities
//--------------------------------------------------------------------------------
void MakeDirectory(const std::string &path)
{
    if (gSystem->AccessPathName(path.c_str())) {
        std::string cmd = "mkdir -p " + path;
        int ret = gSystem->Exec(cmd.c_str());
        if (ret != 0) {
            std::cerr << "WARNING: Could not create directory: " << path
            << " (error code " << ret << ")" << std::endl;
        }
    }
}
// Mark known-bad boards
bool isBadBoard(int sector, int ib)
{
  return ((sector==50 && ib==1) ||
          (sector==4  && ib==1)  ||
          (sector==25 && ib==2));
}

// ---------------------------------------------------------------------
// readRunCSV
// ---------------------------------------------------------------------
std::vector<RunInfo> readRunCSV(const std::string &csvPath)
{
  std::vector<RunInfo> out;
  std::ifstream fin(csvPath);
  if(!fin.is_open()){
    std::cerr<<"[ERROR] cannot open runCSV "<< csvPath <<"\n";
    return out;
  }
  // skip header
  {
    std::string hd;
    if(std::getline(fin, hd)) {
      // ignoring
    }
  }
  int lineNum=0;
  while(true){
    std::string line;
    if(!std::getline(fin, line)) break;
    lineNum++;
    if(line.empty()) continue;
    std::stringstream ss(line);
    std::vector<std::string> tokens;
    std::string tk;
    while(std::getline(ss, tk, ',')){
      tokens.push_back(tk);
    }
    if(tokens.size()<3){
      std::cerr<<"[WARN] line#"<< lineNum <<" => not 3 columns => skip\n";
      continue;
    }
    RunInfo ri;
    ri.runNumber    = std::stoi(tokens[0]);
    ri.width        = std::stoi(tokens[1]);
    std::string off = tokens[2];   // "original" or numeric
    if(off=="original"){
      ri.isOriginal    = true;
      ri.numericOffset = std::nan("");
    } else {
      ri.isOriginal    = false;
      try {
        ri.numericOffset = std::stod(off);
      } catch(...){
        ri.numericOffset = std::nan("");
        ri.isOriginal    = true;
      }
    }
    out.push_back(ri);
  }
  fin.close();
  return out;
}

// ---------------------------------------------------------------------
// readVopCSV => parse tower defaults
// ---------------------------------------------------------------------
std::vector<TowerDefault> readVopCSV(const std::string &csvPath)
{
  std::vector<TowerDefault> out;
  std::ifstream fin(csvPath);
  if(!fin.is_open()){
    std::cerr<<"[ERROR] cannot open vopCSV="<< csvPath <<"\n";
    return out;
  }
  // skip header
  {
    std::string hd;
    if(std::getline(fin, hd)){
      // ignoring
    }
  }
  int lineNo=0;
  while(true){
    std::string line;
    if(!std::getline(fin, line)) break;
    lineNo++;
    if(line.empty()) continue;
    std::stringstream ss(line);
    std::vector<std::string> tokens;
    std::string tk;
    while(std::getline(ss, tk, ',')){
      tokens.push_back(tk);
    }
    // Expect at least 8 columns
    // #1 => sector, #2=>ib, #4=>iphi, #5=>ieta, #7=> defaultOffset
    if(tokens.size()<8){
      std::cerr<<"[WARN] line #"<< lineNo <<" => <8 columns => skip\n";
      continue;
    }
    TowerDefault td;
    td.sector       = std::stoi(tokens[1]);
    td.ib           = std::stoi(tokens[2]);
    td.iphi         = std::stoi(tokens[4]);
    td.ieta         = std::stoi(tokens[5]);
    td.defaultOffset= 0.;
    try {
      td.defaultOffset= std::stod(tokens[7]);
    } catch(...){
      td.defaultOffset= 0.;
    }
    out.push_back(td);
  }
  fin.close();
  return out;
}

void zeroOutBadBoards(TH2F* h2)
{
  if(!h2) return;
  // if bad => set bin content = -9999
  for(int s=0; s<64; s++){
    for(int b=0; b<6; b++){
      if(!isBadBoard(s,b)) continue;
      int localSec= s%32;
      int phiMin= localSec*8;
      int phiMax= phiMin+7;
      int etaMin=0, etaMax=0;
      if(s<32){
        etaMin=48 + 8*b;
        etaMax= etaMin+7;
      } else {
        int base= 40 - 8*b;
        etaMin= base;
        etaMax= base+7;
      }
      for(int xx=phiMin; xx<=phiMax; xx++){
        for(int yy=etaMin; yy<=etaMax; yy++){
          int ibx= h2->GetXaxis()->FindBin(xx);
          int iby= h2->GetYaxis()->FindBin(yy);
          h2->SetBinContent(ibx,iby,-9999.);
        }
      }
    }
  }
}

// ---------------------------------------------------------------------
// getH2CEMC => open test-<run>.root, read "h2CEMC"
// ---------------------------------------------------------------------
TH2F* getH2CEMC(const std::string &dir, int runID)
{
  char fname[512];
  std::snprintf(fname,sizeof(fname),"%s/test-%d.root", dir.c_str(), runID);
  TFile *f= TFile::Open(fname,"READ");
  if(!f || f->IsZombie()){
    if(f) {f->Close(); delete f;}
    std::cerr<<"[WARN] cannot open "<<fname<<"\n";
    return nullptr;
  }
  TH2F* h2= dynamic_cast<TH2F*>( f->Get("h2CEMC") );
  if(!h2){
    std::cerr<<"[WARN] no h2CEMC in "<<fname<<"\n";
    f->Close(); delete f;
    return nullptr;
  }
  TH2F* out= (TH2F*) h2->Clone(Form("h2CEMC_%d", runID));
  out->SetDirectory(0);
  f->Close(); delete f;
  return out;
}


void buildBestPulseWidthMap(const std::string &oldRunCSV,
                            const std::string &oldInputDir,
                            const std::string &newRunCSV,
                            const std::string &newInputDir,
                            const std::string &vopCSV,
                            // OUT references => your final maps
                            std::map< std::pair<int,int>, BoardFitData > &oldFitMap,
                            std::map< std::pair<int,int>, BoardFitData > &newFitMap)
{
    std::cout << "\n\033[1;92m==============================================\033[0m" << std::endl;
    std::cout <<   "\033[1;92m=== Starting buildBestPulseWidthMap(...)  ===\033[0m" << std::endl;
    std::cout <<   "\033[1;92m==============================================\033[0m\n" << std::endl;

    // =====================================================================
    // 1) read old/new CSV => build runDB
    // =====================================================================
    std::cout << "\033[1;34m[STEP 1] Reading old/new CSV files and building runDB...\033[0m" << std::endl;

    std::vector<RunInfo> oldRuns = readRunCSV(oldRunCSV);
    std::vector<RunInfo> newRuns = readRunCSV(newRunCSV);

    std::map<int, RunInfo> runDB;
    std::set<int> oldSet;

    for (auto &rr : oldRuns) {
        runDB[ rr.runNumber ] = rr;
        oldSet.insert(rr.runNumber);
    }
    for (auto &rr : newRuns) {
        runDB[ rr.runNumber ] = rr;
    }
    auto isOld = [&](int runID){ return (oldSet.count(runID) > 0); };

    std::cout << "\033[1;36m   -> oldRuns count: " << oldRuns.size()
              << ", newRuns count: " << newRuns.size() << "\033[0m" << std::endl;


    // =====================================================================
    // 2) read VOP => build (sector,ib) map
    // =====================================================================
    std::cout << "\033[1;34m[STEP 2] Reading VOP CSV and building (sector, ib) map...\033[0m" << std::endl;

    std::vector<TowerDefault> allTows = readVopCSV(vopCSV);
    std::map<std::pair<int,int>, bool> sectorIBmap;
    for(auto &td : allTows){
        sectorIBmap[{td.sector, td.ib}] = true;
    }
    std::cout << "\033[1;36m   -> Total towers read: " << allTows.size()
              << ", Unique (sector, ib) pairs: " << sectorIBmap.size() << "\033[0m" << std::endl;


    // =====================================================================
    // 3) We'll store data in dataMap[0 or 1][(sector,ib)][width]
    //    idx=0 => old, idx=1 => new
    // =====================================================================
    std::cout << "\033[1;34m[STEP 3] Preparing data structures to store run data...\033[0m" << std::endl;

    std::map<int, std::map< std::pair<int,int>, std::map<int, std::vector<XYerr>> >> dataMap;

    // Helper => fill data from a single run
    std::cout << "\033[1;34m[INFO] Defining helper lambda for processing runs...\033[0m" << std::endl;
    auto processRun = [&](int runID, const std::string &dir, bool oldNotNew){
        auto it = runDB.find(runID);
        if(it==runDB.end()) {
            std::cout << "\033[1;33m   -> WARNING: runID " << runID
                      << " not found in runDB, skipping.\033[0m" << std::endl;
            return;
        }

        RunInfo &ri = it->second;
        TH2F* h2 = getH2CEMC(dir, runID);
        if(!h2) {
            std::cout << "\033[1;33m   -> WARNING: TH2F not found for runID "
                      << runID << ", skipping.\033[0m" << std::endl;
            return;
        }
        zeroOutBadBoards(h2);

        // for each known board => gather average ADC
        for(auto &sb: sectorIBmap) {
            int s= sb.first.first;
            int b= sb.first.second;
            if(isBadBoard(s,b)) continue;

            // find default offset
            double defOff= 0.;
            {
                double sum=0.;
                int ct=0;
                for(auto &td : allTows){
                    if(td.sector==s && td.ib==b){
                        sum += td.defaultOffset;
                        ct++;
                    }
                }
                defOff = (ct>0? sum/ct : 0.);
            }

            double delta=0.;
            if(!ri.isOriginal){
                delta = ri.numericOffset - defOff;
            }
            // else delta=0 => means "original offset"

            // gather average ADC among the towers in that board
            double sumVal=0., sumErr2=0.;
            int count=0;
            for(auto &td: allTows){
                if(td.sector==s && td.ib==b) {
                    int bx= h2->GetXaxis()->FindBin(td.iphi);
                    int by= h2->GetYaxis()->FindBin(td.ieta);
                    double val = h2->GetBinContent(bx,by);
                    double err = h2->GetBinError(bx,by);
                    if(val<0) continue; // skip invalid
                    sumVal += val;
                    sumErr2+=(err*err);
                    count++;
                }
            }
            if(count<1) continue;

            double meanADC = sumVal / count;
            double meanErr = std::sqrt(sumErr2)/ count;
            if(meanADC <= 0) continue;

            // store it
            XYerr ex{ delta, meanADC, meanErr };
            int idx= (oldNotNew? 0 : 1);

            dataMap[idx][{s,b}][ ri.width ].push_back(ex);
        }

        delete h2;
    };

    // =====================================================================
    // 3A) Process old runs
    // =====================================================================
    std::cout << "\033[1;34m[STEP 3A] Processing OLD runs …\033[0m"
              << "  (" << oldRuns.size() << " total)\n";

    size_t idxOld = 0;
    for (const auto &rr : oldRuns)
    {
        ++idxOld;
        std::cout << "   • [OLD " << std::setw(3) << idxOld << "/"
                  << oldRuns.size() << "]  run=" << rr.runNumber << "  … ";
        std::cout.flush();

        processRun(rr.runNumber, oldInputDir, /*oldNotNew=*/true);

        std::cout << "\033[32mDONE\033[0m\n";
    }
    
    // 3B) Process new runs
    std::cout << "\033[1;34m[STEP 3B] Processing NEW runs …\033[0m"
              << "  (" << newRuns.size() << " total)\n";

    size_t idxNew = 0;
    for (const auto &rr : newRuns)
    {
        ++idxNew;
        std::cout << "   • [NEW " << std::setw(3) << idxNew << "/"
                  << newRuns.size() << "]  run=" << rr.runNumber << "  … ";
        std::cout.flush();

        processRun(rr.runNumber, newInputDir, /*oldNotNew=*/false);

        std::cout << "\033[32mDONE\033[0m\n";
    }


    // =====================================================================
    // 4) find bestW for OLD => offset=0 near 3000
    // =====================================================================
    std::cout << "\033[1;34m[STEP 4] Finding best pulse width for OLD boards (offset=0 near 3000)...\033[0m" << std::endl;

    std::map<std::pair<int,int>,int> bestWmap;
    {
        // offset=0 => store ADC
        std::map< std::pair<int,int>, std::map<int,std::vector<double>> > zeroStore;
        for(auto &sbItem : dataMap[0]) // idx=0 => old
        {
            auto &Wmap = sbItem.second; // Wmap => { width => vector<XYerr> }
            for(auto &wItem : Wmap){
                int w = wItem.first;
                for(auto &xy: wItem.second){
                    // only if x ~ 0 => offset=0
                    if(std::fabs(xy.x) < 1e-9 && xy.y>0){
                        zeroStore[ sbItem.first ][ w ].push_back(xy.y);
                    }
                }
            }
        }
        // now pick which width's mean ADC is closest to 3000
        for(auto &zz : zeroStore){
            int s= zz.first.first;
            int b= zz.first.second;
            double bestDiff=1e9;
            int bestW=-1;
            for(auto &ww: zz.second){
                int w= ww.first;
                auto &vals = ww.second;
                if(vals.empty()) continue;
                double sum=0.;
                for(auto v: vals) sum+= v;
                double mean= sum / vals.size();
                double df= std::fabs(mean - 3000.);
                if(df < bestDiff){
                    bestDiff= df;
                    bestW= w;
                }
            }
            bestWmap[{s,b}] = bestW;
        }
    }
    std::cout << "\033[1;36m   -> Computed bestW for "
              << bestWmap.size() << " (sector, board) combos (OLD)\033[0m" << std::endl;


    // =====================================================================
    // 5) find candidate widths for NEW => offset=0 near 3000, but we store
    //    *all* widths sorted by closeness. We'll pick the first valid slope
    //    in step6.
    // =====================================================================
    std::cout << "\033[1;34m[STEP 5] Finding candidate pulse widths for NEW boards (offset=0 near 3000)...\033[0m" << std::endl;

    // We'll store a sorted list of candidate widths for each (sector,board),
    // from "closest to 3000" to "furthest".
    std::map< std::pair<int,int>, std::vector<int> > candidatesNewBW;
    {
        // offset=0 => gather ADC values
        std::map< std::pair<int,int>, std::vector< std::pair<int,double> >> rawNew;
        // rawNew[(s,b)] => vector of (width,meanADC) for offset=0

        for(auto &sbN: dataMap[1]) // idx=1 => new
        {
            auto &Wmap= sbN.second;
            for(auto &kv: Wmap){
                int w= kv.first;
                // gather all points with offset=0
                double sum=0.;
                int nPts=0;
                for(auto &xy : kv.second){
                    if(std::fabs(xy.x)<1e-9 && xy.y>0){
                        sum += xy.y;
                        nPts++;
                    }
                }
                if(nPts>0){
                    double meanADC= sum/nPts;
                    rawNew[sbN.first].push_back({ w, meanADC });
                }
            }
        }

        // now sort them by closeness to 3000
        for(auto &xx : rawNew){
            int s= xx.first.first;
            int b= xx.first.second;

            std::sort(xx.second.begin(), xx.second.end(),
                      [&](const auto &a, const auto &b){
                          double diffA = std::fabs(a.second - 3000.);
                          double diffB = std::fabs(b.second - 3000.);
                          return diffA < diffB;
                      });

            for(auto &wm : xx.second){
                candidatesNewBW[{s,b}].push_back(wm.first);
            }
        }
    }
    std::cout << "\033[1;36m   -> Prepared sorted candidate widths for NEW boards: "
              << candidatesNewBW.size() << "\033[0m" << std::endl;


    // =====================================================================
    // 6) Exponential fit => store in oldFitMap / newFitMap
    //    We keep the old approach for OLD boards (single bestW),
    //    but for NEW boards, we attempt each candidate in ascending distance
    //    from 3000. If the slope is negative, skip it & log the attempt.
    // =====================================================================
    std::cout << "\033[1;34m[STEP 6] Fitting data with exponential model and filling result maps...\033[0m" << std::endl;

    // We'll store a more detailed record so each negative attempt includes
    // the final "bestW" (if eventually found) and final slope.
    struct FullNegReport {
        int s, ib;
        int negWidth;
        double negSlope;
        int finalW;      // the bestW that was ultimately chosen
        double finalK;   // final slope for that bestW
    };
    std::vector<FullNegReport> negSlopesNewDetail;

    // A sub-function to do the actual fit
    auto doExponentialFit = [&](bool oldNotNew,
                                int sector, int ib,
                                int bestW) -> BoardFitData
    {
        BoardFitData result;
        result.isOld   = oldNotNew;
        result.sector  = sector;
        result.ib      = ib;
        result.bestW   = bestW;

        int idx = (oldNotNew ? 0 : 1);
        auto &arr = dataMap[idx][{sector, ib}][bestW];
        if(arr.empty()) {
            // no data => empty BoardFitData
            return result;
        }

        // copy all widths => so we have them in .data
        for(auto &kv: dataMap[idx][{sector, ib}]){
            result.data[kv.first] = kv.second;
        }

        // build TGraphErrors for the bestW's data points
        TGraphErrors *gE = new TGraphErrors((int)arr.size());
        for(int i=0; i<(int)arr.size(); i++){
            gE->SetPoint(i, arr[i].x, arr[i].y);
            gE->SetPointError(i, 0., arr[i].e);
        }

        TF1 fExp("fExp","[0]*exp([1]*x)", -1500,1500);

        // initial guess
        double x0,y0;
        gE->GetPoint(0, x0, y0);
        double yFirst= (y0>0.? y0 : 1000.);
        fExp.SetParameters(yFirst, 0.0005);

        TFitResultPtr fr = gE->Fit(&fExp, "S E R Q");

        double A  = fExp.GetParameter(0);
        double eA = fExp.GetParError(0);
        double k  = fExp.GetParameter(1);
        double eK = fExp.GetParError(1);

        result.amplitude     = A;
        result.amplitudeErr  = eA;
        result.kParam        = k;
        result.kParamErr     = eK;

        delete gE;
        return result;
    };


    // --- old => fill oldFitMap
    int oldFitCount = 0;
    std::map<int,int> oldWidthUsage;

    for (auto &pr : dataMap[0]) {
        int s= pr.first.first;
        int b= pr.first.second;
        int w= bestWmap[{s,b}];

        // track widths usage
        for (auto &kv : pr.second) {
            int thisW = kv.first;
            oldWidthUsage[thisW] += 1;
        }

        BoardFitData bfd = doExponentialFit(/*oldNotNew=*/true, s,b,w);
        oldFitMap[{s,b}] = bfd;
        if(!bfd.data.empty()) {
            oldFitCount++;
        }
    }

    // --- new => fill newFitMap
    int newFitCount = 0;
    std::map<int,int> newWidthUsage;

    for (auto &pr : dataMap[1]) {
        int s= pr.first.first;
        int b= pr.first.second;

        // gather widths usage
        for (auto &kv : pr.second) {
            int thisW = kv.first;
            newWidthUsage[thisW] += 1;
        }

        // We'll pick from candidatesNewBW in ascending closeness to 3000
        auto candIt = candidatesNewBW.find({s,b});
        if(candIt==candidatesNewBW.end() || candIt->second.empty()){
            // no candidate => store empty BoardFitData
            newFitMap[{s,b}] = BoardFitData();
            continue;
        }

        // Attempt each candidate in turn. The first that yields slope>=0 => final
        BoardFitData finalBFD;
        bool foundPositiveSlope=false;

        // We'll store the negative attempts in a local vector, so that once
        // we do or don't find a final slope, we can fill "finalW" in them.
        std::vector<FullNegReport> localNegs;

        for(int bestCandidate : candIt->second){
            BoardFitData testFit = doExponentialFit(/*oldNotNew=*/false, s,b,bestCandidate);

            // if no data => skip
            if(testFit.data.empty())
                continue;

            double slopeK = testFit.kParam;
            if(slopeK >= 0.) {
                // success => adopt
                finalBFD = testFit;
                foundPositiveSlope = true;
                break;
            }
            else {
                // negative => store a partial record for this attempt
                FullNegReport tmp;
                tmp.s = s;
                tmp.ib= b;
                tmp.negWidth= bestCandidate;
                tmp.negSlope= slopeK;
                tmp.finalW  = -1;  // unknown yet
                tmp.finalK  = 0.;
                localNegs.push_back(tmp);
            }
        }

        // Now fill final best
        newFitMap[{s,b}] = finalBFD;
        if(!finalBFD.data.empty()) {
            newFitCount++;
        }

        // If we found a final best => we fill in finalW, finalSlope in all
        // localNegs for this board
        int finalUsedW  = (foundPositiveSlope? finalBFD.bestW  : -1);
        double finalUsedK= (foundPositiveSlope? finalBFD.kParam : 0.);

        for (auto &neg : localNegs) {
            neg.finalW = finalUsedW;
            neg.finalK = finalUsedK;
            negSlopesNewDetail.push_back(neg);
        }
    }

    std::cout << "\033[1;36m   -> Exponential fits completed for OLD boards: " << oldFitCount
              << ", NEW boards: " << newFitCount << "\033[0m" << std::endl;

    // =====================================================================
    // 6B) Print the negative-slope attempts table with final chosen W
    // =====================================================================
    if(!negSlopesNewDetail.empty()) {
        std::cout << "\n\033[1;33m[WARNING] Negative slope encountered for some NEW boards:\033[0m\n";
        // We expand the columns to also show the final chosen width & slope
        std::cout << std::left
                  << std::setw(8) << "Sector"
                  << std::setw(6) << "IB"
                  << std::setw(8) << "negW"
                  << std::setw(12)<< "negSlope"
                  << std::setw(8) << "finalW"
                  << std::setw(12)<< "finalSlope"
                  << "\n";

        // Print a dashed separator
        std::cout << std::string(8+6+8+12+8+12, '-') << "\n";

        for(const auto &ev : negSlopesNewDetail){
            std::cout << std::left
                      << std::setw(8)  << ev.s
                      << std::setw(6)  << ev.ib
                      << std::setw(8)  << ev.negWidth
                      << std::setw(12) << Form("%.4g", ev.negSlope)
                      << std::setw(8)  << ev.finalW
                      << std::setw(12) << Form("%.4g", ev.finalK)
                      << "\n";
        }
        std::cout << std::endl;
    }

    // =====================================================================
    // Final Summary/Check
    // =====================================================================
    std::cout << "\n\033[1;92m=======================================\033[0m" << std::endl;
    std::cout <<   "\033[1;92m===  Summary of buildBestPulseWidthMap ===\033[0m" << std::endl;
    std::cout <<   "\033[1;92m=======================================\033[0m\n" << std::endl;

    // Example tabulated info:
    std::cout << std::left
              << std::setw(25) << " # Old Runs"
              << std::setw(25) << " # New Runs"
              << std::setw(25) << " # Old Fit Boards"
              << std::setw(25) << " # New Fit Boards"
              << std::endl;

    std::cout << std::left
              << std::setw(25) << oldRuns.size()
              << std::setw(25) << newRuns.size()
              << std::setw(25) << oldFitCount
              << std::setw(25) << newFitCount
              << std::endl;

    // ---------------------------------------------------------------------
    // Additional table: how many IBs have each width in data, how many use
    // it as best. Then maybe a short amplitude/k summary if you desire.
    // ---------------------------------------------------------------------
    // 1) find the union of widths observed in oldWidthUsage / newWidthUsage
    std::set<int> allWidths;
    for (auto &kv : oldWidthUsage) allWidths.insert(kv.first);
    for (auto &kv : newWidthUsage) allWidths.insert(kv.first);

    // 2) build counters for “best usage” => oldBestCount[w], newBestCount[w]
    std::map<int,int> oldBestCount, newBestCount;
    for (auto &pr : oldFitMap) {
        int wUsed = pr.second.bestW;
        if (wUsed != 0) {
            oldBestCount[wUsed]++;
        }
    }
    for (auto &pr : newFitMap) {
        int wUsed = pr.second.bestW;
        if (wUsed != 0) {
            newBestCount[wUsed]++;
        }
    }

    std::cout << "\n\033[1;94m=== DETAILED PULSE WIDTH USAGE TABLE (OLD vs NEW) ===\033[0m\n";
    std::cout << std::left
              << std::setw(12) << "Width"
              << std::setw(18) << "OLD(HasData)"
              << std::setw(18) << "OLD(asBEST)"
              << std::setw(18) << "NEW(HasData)"
              << std::setw(18) << "NEW(asBEST)"
              << std::endl;

    // print line
    std::cout << std::string(12+18+18+18+18,'-') << "\n";

    for (auto w : allWidths) {
        int oHas   = (oldWidthUsage.count(w) ? oldWidthUsage[w] : 0);
        int oBest  = (oldBestCount.count(w)  ? oldBestCount[w]  : 0);
        int nHas   = (newWidthUsage.count(w) ? newWidthUsage[w] : 0);
        int nBest  = (newBestCount.count(w)  ? newBestCount[w]  : 0);

        std::cout << std::left
                  << std::setw(12) << w
                  << std::setw(18) << oHas
                  << std::setw(18) << oBest
                  << std::setw(18) << nHas
                  << std::setw(18) << nBest
                  << std::endl;
    }

    // ---------------------------------------------------------------------
    // 3) Optional: list a few example boards showing amplitude / kParam
    //    for old vs new. (Just to confirm everything is good.)
    // ---------------------------------------------------------------------
    std::cout << "\n\033[1;94m=== SAMPLE FIT OUTPUTS (up to 5 boards) ===\033[0m\n";
    int countPrinted=0;
    for (auto &pr : oldFitMap) {
        if (countPrinted>=5) break;
        int s= pr.first.first;
        int b= pr.first.second;
        const BoardFitData &bfdO = pr.second;

        // see if newFitMap has that board
        auto itN = newFitMap.find({s,b});
        if (itN == newFitMap.end()) continue;

        const BoardFitData &bfdN = itN->second;
        if (bfdO.data.empty() || bfdN.data.empty()) continue;

        std::cout << "  Board => s="<< s <<", ib="<< b
                  << " | Old: A="<<bfdO.amplitude<<", k="<<bfdO.kParam
                  << " || New: A="<<bfdN.amplitude<<", k="<<bfdN.kParam
                  << std::endl;
        countPrinted++;
    }


    std::cout << "\n\033[1;92m=== buildBestPulseWidthMap COMPLETED SUCCESSFULLY ===\033[0m\n" << std::endl;


    // --------------------------------------------------------------------------------
    // ADDITIONAL STEP: Produce two CSV files with full details for OLD and NEW data.
    // --------------------------------------------------------------------------------
    const std::string outDir   = "/Users/patsfan753/Desktop/PMTgainsAna/output";
    gSystem->Exec(Form("mkdir -p %s", outDir.c_str()));

    std::string oldCSVpath = outDir + "/oldPulseWidthData.csv";
    std::string newCSVpath = outDir + "/newPulseWidthData.csv";

    std::ofstream oldCSV(oldCSVpath);
    std::ofstream newCSV(newCSVpath);

    if(!oldCSV.is_open() || !newCSV.is_open()){
        std::cerr << "\033[1;31m[ERROR] Could not open the CSV files for writing:\033[0m\n"
                  << "  " << oldCSVpath << "\n"
                  << "  " << newCSVpath << "\n";
    }
    else {
        // Write headers
        oldCSV << "sector,ib,bestW,amplitude,amplitudeErr,kParam,kParamErr,width,offset,avgADC,adcErr\n";
        newCSV << "sector,ib,bestW,amplitude,amplitudeErr,kParam,kParamErr,width,offset,avgADC,adcErr\n";

        // Fill OLD CSV
        for(const auto &pr : oldFitMap){
            int s= pr.first.first;
            int b= pr.first.second;
            const BoardFitData &bfd = pr.second;
            // for each width => for each XYerr
            for(const auto &kv : bfd.data){
                int w = kv.first;
                const auto &vecErrs = kv.second;
                for(const auto &xy : vecErrs){
                    oldCSV
                      << s << ","
                      << b << ","
                      << bfd.bestW << ","
                      << bfd.amplitude << ","
                      << bfd.amplitudeErr << ","
                      << bfd.kParam << ","
                      << bfd.kParamErr << ","
                      << w << ","
                      << xy.x << ","
                      << xy.y << ","
                      << xy.e << "\n";
                }
            }
        }

        // Fill NEW CSV
        for(const auto &pr : newFitMap){
            int s= pr.first.first;
            int b= pr.first.second;
            const BoardFitData &bfd = pr.second;
            if(bfd.data.empty()) continue;  // skip empty fits
            for(const auto &kv : bfd.data){
                int w = kv.first;
                const auto &vecErrs = kv.second;
                for(const auto &xy : vecErrs){
                    newCSV
                      << s << ","
                      << b << ","
                      << bfd.bestW << ","
                      << bfd.amplitude << ","
                      << bfd.amplitudeErr << ","
                      << bfd.kParam << ","
                      << bfd.kParamErr << ","
                      << w << ","
                      << xy.x << ","
                      << xy.y << ","
                      << xy.e << "\n";
                }
            }
        }

        oldCSV.close();
        newCSV.close();

        std::cout << "\033[1;92m[INFO] Successfully wrote CSV summaries:\033[0m\n"
                  << "  " << oldCSVpath << "\n"
                  << "  " << newCSVpath << "\n";
    }
}



void produceCompare32vsBest(
    const std::map< std::pair<int,int>, BoardFitData > &oldFitMap,
    const std::map< std::pair<int,int>, BoardFitData > &newFitMap
)
{
    std::cout << "\n[DEBUG] produceCompare32vsBest() => building hist: old(32,0) vs old(best), new(32,0) vs new(best)\n";

    //----------------------------------------------------------------------
    // 1) Create four histograms:
    //      - histOld32, histOldBest
    //      - histNew32, histNewBest
    //----------------------------------------------------------------------
    TH1F* histOld32   = new TH1F("histOld32",   "OLD=width32, offset=0;ADC;Count",128,0,16384);
    histOld32->SetLineColor(kBlue);
    histOld32->SetLineWidth(3);

    TH1F* histOldBest = new TH1F("histOldBest","OLD bestW offset=0;ADC;Count",128,0,16384);
    histOldBest->SetLineColor(kMagenta+2);
    histOldBest->SetLineWidth(3);

    TH1F* histNew32   = new TH1F("histNew32",   "NEW=width32, offset=0;ADC;Count",128,0,16384);
    histNew32->SetLineColor(kBlue);
    histNew32->SetLineWidth(3);

    TH1F* histNewBest = new TH1F("histNewBest", "NEW bestW offset=0;ADC;Count",128,0,16384);
    histNewBest->SetLineColor(kRed);
    histNewBest->SetLineWidth(3);

    // We also track which boards contribute to each ADC bin for "new bestW"
    // so we can print out boards that cause spikes in [2500..3500].
    std::map<int, std::vector< std::pair<int,int> >> newBestContributors;

    //----------------------------------------------------------------------
    // 2) Fill OLD sample histograms => width=32 offset=0, and bestW offset=0
    //----------------------------------------------------------------------
    for (auto &pr : oldFitMap)
    {
        const BoardFitData &bfd = pr.second;

        // (A) bestW => fill histOldBest
        int wBest = bfd.bestW;  // already computed in buildBestPulseWidthMap(...)
        if (wBest > 0) {
            auto itVec = bfd.data.find(wBest);
            if (itVec != bfd.data.end()) {
                for (auto &xy : itVec->second) {
                    // We want only offset=0 => x near 0
                    if (std::fabs(xy.x) < 1e-9 && xy.y > 0.) {
                        histOldBest->Fill(xy.y);
                    }
                }
            }
        }

        // (B) width=32 => fill histOld32
        auto it32 = bfd.data.find(32);
        if (it32 != bfd.data.end()) {
            for (auto &xy : it32->second) {
                if (std::fabs(xy.x) < 1e-9 && xy.y > 0.) {
                    histOld32->Fill(xy.y);
                }
            }
        }
    }

    //----------------------------------------------------------------------
    // 3) Fill NEW sample histograms => width=32 offset=0, and bestW offset=0
    //----------------------------------------------------------------------
    for (auto &pr : newFitMap)
    {
        const BoardFitData &bfd = pr.second;
        int s      = bfd.sector;
        int b      = bfd.ib;
        int wBest  = bfd.bestW;  // from buildBestPulseWidthMap(...)

        // (A) bestW => fill histNewBest
        if (wBest > 0) {
            auto itVec = bfd.data.find(wBest);
            if (itVec != bfd.data.end()) {
                for (auto &xy : itVec->second) {
                    if (std::fabs(xy.x) < 1e-9 && xy.y>0.) {
                        double adcVal = xy.y;
                        histNewBest->Fill(adcVal);

                        // record => which boards contributed to this bin
                        int binIndex = histNewBest->FindBin(adcVal);
                        newBestContributors[binIndex].push_back({ s, b });
                    }
                }
            }
        }

        // (B) width=32 => fill histNew32
        auto it32 = bfd.data.find(32);
        if (it32 != bfd.data.end()) {
            for (auto &xy : it32->second) {
                if (std::fabs(xy.x) < 1e-9 && xy.y>0.) {
                    histNew32->Fill(xy.y);
                }
            }
        }
    }

    //----------------------------------------------------------------------
    // 4) Print out boards that cause "spikes" in NEW bestW, ADC in [2500..3500]
    //----------------------------------------------------------------------
    std::cout << "\n[DEBUG] => Listing boards that cause spikes in NEW bestW, ADC in [2500..3500]...\n";
    {
        int firstBin = histNewBest->FindBin(2500.);
        int lastBin  = histNewBest->FindBin(3500.);
        for(int ibin = firstBin; ibin <= lastBin; ibin++){
            int count = (int) histNewBest->GetBinContent(ibin);
            if (count < 1) continue;

            double center = histNewBest->GetBinCenter(ibin);
            std::cout << "  Bin #" << ibin << " => center=" << center
                      << ", count=" << count << ":\n";

            auto &boards = newBestContributors[ibin];
            for (auto &sb : boards) {
                int sec = sb.first;
                int ibd = sb.second;
                std::cout << "     sector=" << sec
                          << ", IB=" << ibd << "\n";
            }
        }
    }

    //----------------------------------------------------------------------
    // 5) A small helper => binStats(...) returns (#below, #inRange, #above, #total)
    //----------------------------------------------------------------------
    auto binStats = [&](TH1F *hist){
        int below=0, inRange=0, above=0;
        int nbins= hist->GetNbinsX();
        for(int i=1; i<=nbins; i++){
            double xC = hist->GetBinCenter(i);
            int    c  = (int) hist->GetBinContent(i);
            if      (xC < 2500.) below   += c;
            else if (xC > 3500.) above   += c;
            else                 inRange += c;
        }
        int total = below + inRange + above;
        return std::make_tuple(below, inRange, above, total);
    };

    //----------------------------------------------------------------------
    // 6) Build final 2×1 TCanvas => "compareADC_32_vs_best.png"
    //    with top sub‐pads (hist + lines) and bottom sub‐pads (ratio).
    //----------------------------------------------------------------------
    TCanvas cF("cF","compareADC_32_vs_best(OLD vs NEW)",1800,700);
    cF.Divide(2,1);

    //---------------------------------------------------------
    // LEFT side => old(32 vs best)
    //---------------------------------------------------------
    cF.cd(1);
    TPad* subPad1= (TPad*)gPad;
    subPad1->Divide(1,2,0.001,0.001);

    // (TOP) => old
    subPad1->cd(1);
    double mRef  = histOld32->GetMaximum();
    double mBest = histOldBest->GetMaximum();
    double yMax  = std::max(mRef, mBest) * 1.2;

    histOld32->SetMaximum(yMax);
    histOld32->SetTitle("OLD Data: width=32 vs. best (#Delta=0)");
    histOld32->GetXaxis()->SetTitle("Average ADC");
    histOld32->GetYaxis()->SetTitle("Count of IBs");
    histOld32->Draw("HIST");
    histOldBest->Draw("HIST SAME");

    // green lines at [2500..3500]
    TLine lA(2500,0,2500,yMax);
    lA.SetLineColor(kGreen+2);
    lA.SetLineStyle(2);
    lA.SetLineWidth(2);
    lA.Draw("same");

    TLine lB(3500,0,3500,yMax);
    lB.SetLineColor(kGreen+2);
    lB.SetLineStyle(2);
    lB.SetLineWidth(2);
    lB.Draw("same");

    // small legend for the top-left
    TLegend legO(0.42,0.70,0.77,0.88);
    legO.AddEntry(histOld32,   "OLD Sample (TP=32)",   "l");
    legO.AddEntry(histOldBest, "OLD Sample (best TP)", "l");
    legO.AddEntry(&lA,         "Range: [2500..3500]",  "l");
    legO.Draw();

    // print tallies for old(32) and old(best)
    {
        auto [b32, i32, a32, t32]        = binStats(histOld32);
        auto [bBest, iBest, aBest, tBest]= binStats(histOldBest);

        std::cout << "\n[DEBUG] OLD(32,0) => total IB="<< t32
                  << " => below="<< b32
                  << ", inRange="<< i32
                  << ", above="<< a32 <<"\n";
        if(t32>0){
            std::cout << "         => % below= "<< (100.*b32/t32)
                      <<", % inRange= "<< (100.*i32/t32)
                      <<", % above= "<< (100.*a32/t32) <<"\n";
        }
        std::cout << "[DEBUG] OLD(best,0) => total IB="<< tBest
                  << " => below="<< bBest
                  << ", inRange="<< iBest
                  << ", above="<< aBest <<"\n";
        if(tBest>0){
            std::cout << "         => % below= "<< (100.*bBest/tBest)
                      <<", % inRange= "<< (100.*iBest/tBest)
                      <<", % above= "<< (100.*aBest/tBest) <<"\n";
        }

        // small text blocks on the left side
        TLatex latL, latR;
        latL.SetNDC(true);
        latR.SetNDC(true);
        latL.SetTextSize(0.04);
        latR.SetTextSize(0.04);

        double yLeft = 0.6;
        latL.SetTextColor(kBlack);
        latL.DrawLatex(0.35, yLeft, Form("OLD(32) total IB=%d", t32));
        yLeft -= 0.05;
        if(t32>0){
            latL.SetTextColor(kBlue);
            latL.DrawLatex(0.35, yLeft, Form("below2500=%.1f%%",100.*b32/t32));
            yLeft -= 0.05;

            latL.SetTextColor(kGreen+2);
            latL.DrawLatex(0.35, yLeft, Form("[2500..3500]=%.1f%%",100.*i32/t32));
            yLeft -= 0.05;

            latL.SetTextColor(kRed);
            latL.DrawLatex(0.35, yLeft, Form("above3500=%.1f%%",100.*a32/t32));
            yLeft -= 0.05;
        }

        double yRight= 0.6;
        latR.SetTextColor(kBlack);
        latR.DrawLatex(0.63, yRight, Form("OLD(best) total IB=%d", tBest));
        yRight -= 0.05;
        if(tBest>0){
            latR.SetTextColor(kBlue);
            latR.DrawLatex(0.63,yRight,Form("below2500=%.1f%%",100.*bBest/tBest));
            yRight -= 0.05;

            latR.SetTextColor(kGreen+2);
            latR.DrawLatex(0.63,yRight,Form("[2500..3500]=%.1f%%",100.*iBest/tBest));
            yRight -= 0.05;

            latR.SetTextColor(kRed);
            latR.DrawLatex(0.63,yRight,Form("above3500=%.1f%%",100.*aBest/tBest));
            yRight -= 0.05;
        }
    }

    // (BOTTOM) => ratio old(best)/old(32)
    subPad1->cd(2);
    TH1F* hRatioO = (TH1F*) histOldBest->Clone("ratioOld");
    hRatioO->Divide(histOld32);
    hRatioO->SetTitle("");
    hRatioO->GetXaxis()->SetTitle("ADC");
    hRatioO->GetYaxis()->SetTitle("OLD(best)/OLD(32)");
    hRatioO->GetYaxis()->SetNdivisions(505);

    double rmO = 0.;
    for(int i=1; i<= hRatioO->GetNbinsX(); i++){
        double val= hRatioO->GetBinContent(i);
        double err= hRatioO->GetBinError(i);
        if (val+err > rmO) rmO = val+err;
    }
    hRatioO->SetMaximum(std::max(3.0, 1.2*rmO));
    hRatioO->SetMinimum(0.);
    hRatioO->Draw("E0");

    //---------------------------------------------------------
    // RIGHT side => new(32 vs best)
    //---------------------------------------------------------
    cF.cd(2);
    TPad* subPad2= (TPad*)gPad;
    subPad2->Divide(1,2,0.001,0.001);

    // (TOP) => new
    subPad2->cd(1);
    double mm1 = histNew32->GetMaximum();
    double mm2 = histNewBest->GetMaximum();
    double yM  = std::max(mm1, mm2)*1.2;
    histNew32->SetMaximum(yM);

    histNew32->SetTitle("NEW Sample: width=32 vs. best (#Delta=0)");
    histNew32->GetXaxis()->SetTitle("Average ADC");
    histNew32->GetYaxis()->SetTitle("Count of IBs");
    histNew32->Draw("HIST");
    histNewBest->Draw("HIST SAME");

    // green lines at [2500..3500]
    TLine LN1(2500,0,2500,yM);
    LN1.SetLineColor(kGreen+2);
    LN1.SetLineStyle(2);
    LN1.SetLineWidth(2);
    LN1.Draw("same");

    TLine LN2(3500,0,3500,yM);
    LN2.SetLineColor(kGreen+2);
    LN2.SetLineStyle(2);
    LN2.SetLineWidth(2);
    LN2.Draw("same");

    TLegend legN(0.42,0.70,0.77,0.88);
    legN.AddEntry(histNew32,   "NEW Sample (TP=32)",   "l");
    legN.AddEntry(histNewBest, "NEW Sample (best TP)", "l");
    legN.AddEntry(&LN1,        "Range: [2500..3500]",  "l");
    legN.Draw();

    // print tallies for new(32) and new(best)
    {
        auto [b32,i32,a32,t32]        = binStats(histNew32);
        auto [bBest,iBest,aBest,tBest]= binStats(histNewBest);

        std::cout << "\n[DEBUG] NEW(32,0) => total IB="<< t32
                  << " => below="<< b32
                  << ", inRange="<< i32
                  << ", above="<< a32 <<"\n";
        if(t32>0){
            std::cout << "         => % below= "<< (100.*b32/t32)
                      <<", % inRange= "<< (100.*i32/t32)
                      <<", % above= "<< (100.*a32/t32) <<"\n";
        }
        std::cout << "[DEBUG] NEW(best,0) => total IB="<< tBest
                  << " => below="<< bBest
                  << ", inRange="<< iBest
                  << ", above="<< aBest <<"\n";
        if(tBest>0){
            std::cout << "         => % below= "<< (100.*bBest/tBest)
                      <<", % inRange= "<< (100.*iBest/tBest)
                      <<", % above= "<< (100.*aBest/tBest) <<"\n";
        }

        TLatex latL, latR;
        latL.SetNDC(true);
        latR.SetNDC(true);
        latL.SetTextSize(0.04);
        latR.SetTextSize(0.04);

        double yLeft= 0.6;
        latL.SetTextColor(kBlack);
        latL.DrawLatex(0.35,yLeft, Form("NEW(32) total IB=%d", t32));
        yLeft -= 0.05;
        if(t32>0){
            latL.SetTextColor(kBlue);
            latL.DrawLatex(0.35,yLeft,Form("below2500=%.1f%%",100.*b32/t32));
            yLeft -= 0.05;

            latL.SetTextColor(kGreen+2);
            latL.DrawLatex(0.35,yLeft,Form("[2500..3500]=%.1f%%",100.*i32/t32));
            yLeft -= 0.05;

            latL.SetTextColor(kRed);
            latL.DrawLatex(0.35,yLeft,Form("above3500=%.1f%%",100.*a32/t32));
            yLeft -= 0.05;
        }

        double yRight= 0.6;
        latR.SetTextColor(kBlack);
        latR.DrawLatex(0.63,yRight, Form("NEW(best) total IB=%d", tBest));
        yRight -= 0.05;
        if(tBest>0){
            latR.SetTextColor(kBlue);
            latR.DrawLatex(0.63,yRight,Form("below2500=%.1f%%",100.*bBest/tBest));
            yRight -= 0.05;

            latR.SetTextColor(kGreen+2);
            latR.DrawLatex(0.63,yRight,Form("[2500..3500]=%.1f%%",100.*iBest/tBest));
            yRight -= 0.05;

            latR.SetTextColor(kRed);
            latR.DrawLatex(0.63,yRight,Form("above3500=%.1f%%",100.*aBest/tBest));
            yRight -= 0.05;
        }
    }

    // (BOTTOM) => ratio new(best)/new(32)
    subPad2->cd(2);
    TH1F* hRatioN= (TH1F*) histNewBest->Clone("ratioNew");
    hRatioN->Divide(histNew32);
    hRatioN->SetTitle("");
    hRatioN->GetXaxis()->SetTitle("ADC");
    hRatioN->GetYaxis()->SetTitle("NEW(best)/NEW(32)");
    hRatioN->GetYaxis()->SetNdivisions(505);

    double rmN=0.;
    for(int i=1; i<= hRatioN->GetNbinsX(); i++){
        double val= hRatioN->GetBinContent(i);
        double err= hRatioN->GetBinError(i);
        if (val+err > rmN) rmN= val+err;
    }
    hRatioN->SetMaximum(std::max(3.0,1.2*rmN));
    hRatioN->SetMinimum(0.);
    hRatioN->Draw("E0");

    //----------------------------------------------------------------------
    // 7) Save final figure
    //----------------------------------------------------------------------
    std::string outTwo = "/Users/patsfan753/Desktop/PMTgainsAna/output/compareADC_32_vs_best.png";
    cF.SaveAs(outTwo.c_str());
    std::cout << "[INFO] => Created " << outTwo << " via produceCompare32vsBest()\n";
}




void produceColorCodedMapsAndExtra(
    const std::map< std::pair<int,int>, BoardFitData > &oldFitMap,
    const std::map< std::pair<int,int>, BoardFitData > &newFitMap,
    const std::string &summary
)
{
    //----------------------------------------------------------------------
    // 0) Collect all boards from oldFitMap + newFitMap
    //----------------------------------------------------------------------
    std::set<std::pair<int,int>> allBoards;
    // Insert boards from old:
    for (auto &kv : oldFitMap) {
        allBoards.insert(kv.first); // (sector, ib)
    }
    // Insert boards from new:
    for (auto &kv : newFitMap) {
        allBoards.insert(kv.first); // (sector, ib)
    }

    //----------------------------------------------------------------------
    // 1) Helpers to retrieve the “offset=0” ADC from W=32 and from bestW
    //----------------------------------------------------------------------
    auto get32offset0 = [&](bool isOld, int sector, int ib) -> double {
        // pick oldFitMap vs newFitMap
        const auto &refMap = (isOld ? oldFitMap : newFitMap);
        auto it = refMap.find({sector, ib});
        if (it == refMap.end()) return -1.0;

        const BoardFitData &bfd = it->second;
        auto itW = bfd.data.find(32);
        if (itW == bfd.data.end()) return -1.0;

        // look for x=0 => offset=0 => pick y>0
        for (auto &xy : itW->second) {
            if (std::fabs(xy.x) < 1e-9 && xy.y > 0.) {
                return xy.y;
            }
        }
        return -1.0;
    };

    auto getBestOffset0 = [&](bool isOld, int sector, int ib) -> double {
        // pick oldFitMap vs newFitMap
        const auto &refMap = (isOld ? oldFitMap : newFitMap);
        auto it = refMap.find({sector, ib});
        if (it == refMap.end()) return -1.0;

        const BoardFitData &bfd = it->second;
        if (bfd.bestW < 1) return -1.0; // invalid bestW

        auto itW = bfd.data.find(bfd.bestW);
        if (itW == bfd.data.end()) return -1.0;

        for (auto &xy : itW->second) {
            if (std::fabs(xy.x) < 1e-9 && xy.y > 0.) {
                return xy.y;
            }
        }
        return -1.0;
    };

    //----------------------------------------------------------------------
    // 2) We'll keep counters for color-coded categories
    //----------------------------------------------------------------------
    int old_below=0, old_inrange=0, old_above=0, old_none=0;
    int new_below=0, new_inrange=0, new_above=0, new_none=0;

    //----------------------------------------------------------------------
    // 3) A helper to get TBox corners for a (sector, ib) tower-based display
    //----------------------------------------------------------------------
    auto getPadCorners = [&](int s, int b) {
        int localSec = (s % 32);
        double phiMin = 8.*localSec;
        double phiMax = phiMin + 8.;
        double etaMin=0., etaMax=0.;

        if (s < 32) {
            etaMin = 48. + 8.*b;
            etaMax = etaMin + 8.;
        } else {
            etaMin = 40. - 8.*b;
            etaMax = etaMin + 8.;
        }
        if (etaMin > etaMax) std::swap(etaMin, etaMax);
        return std::make_tuple(phiMin, phiMax, etaMin, etaMax);
    };

    //----------------------------------------------------------------------
    // 4) chooseColor(...) => pick color based on ADC, update counters
    //----------------------------------------------------------------------
    auto chooseColor = [&](double adcVal, bool isOld) -> Color_t {
        if (adcVal < 0.) {
            // means no data or invalid => gray
            if (isOld) old_none++; else new_none++;
            return (Color_t)(kGray+1);
        }
        else if (adcVal < 2500.) {
            if (isOld) old_below++; else new_below++;
            return kBlue;
        }
        else if (adcVal > 3500.) {
            if (isOld) old_above++; else new_above++;
            return kRed;
        }
        else {
            if (isOld) old_inrange++; else new_inrange++;
            return (Color_t)(kGreen+2);
        }
    };

    //----------------------------------------------------------------------
    // 5) drawColorMapInPad(...) => draws W=32 or bestW in the *current* pad
    //----------------------------------------------------------------------
    auto drawColorMapInPad = [&](bool isOld, bool w32OrBest) {
        // Iterate over all boards in allBoards set
        for (auto &sb : allBoards) {
            int s = sb.first;
            int b = sb.second;

            double valADC = (w32OrBest
                             ? get32offset0(isOld, s, b)
                             : getBestOffset0(isOld, s, b));
            // corners
            auto [x1, x2, y1, y2] = getPadCorners(s,b);

            TBox *box = new TBox(x1, y1, x2, y2);
            box->SetLineColor(kBlack);
            box->SetLineWidth(1);

            Color_t col = chooseColor(valADC, isOld);
            box->SetFillColor(col);
            box->Draw();

            // label => "s.ib"
            double cx = 0.5*(x1 + x2);
            double cy = 0.5*(y1 + y2);
            TLatex lat;
            lat.SetTextAlign(22);
            lat.SetTextSize(0.017);
            lat.SetTextColor(kBlack);
            lat.DrawLatex(cx, cy, Form("%d.%d", s, b));
        }
    };

    //----------------------------------------------------------------------
    // 6) produceCombinedMap(...) => builds the 2×1 color-coded map
    //    for either old or new => top = W=32, bottom = bestW
    //----------------------------------------------------------------------
    auto produceCombinedMap = [&](bool isOld, const std::string &outFile) {
        std::string cName  = (isOld ? "cMapOldCombined" : "cMapNewCombined");
        std::string cTitle = (isOld
                              ? "Old color map (W=32 on top, bestW bottom)"
                              : "New color map (W=32 on top, bestW bottom)");
        TCanvas cC(cName.c_str(), cTitle.c_str(), 3000, 2400);
        cC.Divide(1, 2, 0.001, 0.001);

        // top => W=32
        cC.cd(1);
        gPad->SetRightMargin(0.20);
        TH2F *hFtop = new TH2F("hFrameTop",";#phi index;#eta index",256,0,256,96,0,96);
        hFtop->SetStats(false);
        hFtop->Draw("COL");
        drawColorMapInPad(isOld, true);

        // build a small legend at the right
        TLegend *legC= new TLegend(0.82, 0.10, 0.98, 0.90);
        legC->SetHeader(isOld
                        ? "OLD: top=W=32, bottom=BEST"
                        : "NEW: top=W=32, bottom=BEST","C");
        legC->SetBorderSize(0);
        legC->SetFillColorAlpha(kWhite,0.6);
        legC->SetTextFont(42);
        legC->SetTextSize(0.035);

        // dummy TBoxes for color legend
        TBox bBlue(0,0,1,1);  bBlue.SetFillColor(kBlue);
        TBox bGreen(0,0,1,1); bGreen.SetFillColor(kGreen+2);
        TBox bRed(0,0,1,1);   bRed.SetFillColor(kRed);
        TBox bGray(0,0,1,1);  bGray.SetFillColor(kGray+1);

        legC->AddEntry(&bBlue, "ADC<2500","f");
        legC->AddEntry(&bGreen,"2500..3500","f");
        legC->AddEntry(&bRed,  "ADC>3500","f");
        legC->AddEntry(&bGray, "No Data","f");
        legC->Draw();

        // bottom => bestW
        cC.cd(2);
        gPad->SetRightMargin(0.20);
        TH2F *hFbot = new TH2F("hFrameBot",";#phi index;#eta index",256,0,256,96,0,96);
        hFbot->SetStats(false);
        hFbot->Draw("COL");
        drawColorMapInPad(isOld, false);

        cC.SaveAs(outFile.c_str());
        std::cout << "[INFO] => Created combined color-coded map: " << outFile << "\n";

        delete hFtop;
        delete hFbot;
    };

    //----------------------------------------------------------------------
    // 7) Actually produce old/new combined color-coded PNG
    //----------------------------------------------------------------------
    std::cout << "[DEBUG] => producing colorMap_old_combined.png & colorMap_new_combined.png\n";
    {
        produceCombinedMap(true, "/Users/patsfan753/Desktop/PMTgainsAna/output/2Dmaps/colorMap_old_combined.png");
        produceCombinedMap(false, "/Users/patsfan753/Desktop/PMTgainsAna/output/2Dmaps/colorMap_new_combined.png");
    }

    //----------------------------------------------------------------------
    // 8) Build the single 2×1 “COLZ” figure for the new sample => W=32 vs bestW
    //----------------------------------------------------------------------
    {
        TCanvas cH("cH","New Sample Heat Maps (COLZ)",2200,1600);
        cH.Divide(1,2);

        // top => new(W=32)
        cH.cd(1);
        gPad->SetRightMargin(0.20);

        TH2F* h2W32 = new TH2F("h2New32","NEW Sample (TPW=32, #Delta=0);#phi index;#eta index",
                               256,0,256,96,0,96);
        h2W32->SetStats(false);
        h2W32->GetZaxis()->SetTitle("avg ADC");

        // fill from get32offset0(false => new)
        for (auto &sb : allBoards) {
            int s= sb.first, b= sb.second;
            double val = get32offset0(false, s, b);
            if (val<0) continue;

            // fill an 8×8 region
            int localSec= (s % 32);
            int phiMin= 8*localSec, phiMax= phiMin+7;
            int etaMin=0, etaMax=0;
            if (s<32) {
                etaMin=48+ 8*b;
                etaMax=etaMin+7;
            } else {
                int base= 40 - 8*b;
                etaMin= base; etaMax= base+7;
            }
            if (etaMin>etaMax) std::swap(etaMin,etaMax);

            for(int xx=phiMin; xx<=phiMax; xx++){
                for(int yy=etaMin; yy<=etaMax; yy++){
                    h2W32->SetBinContent(xx+1, yy+1, val);
                }
            }
        }
        h2W32->Draw("COLZ");

        // also label each IB => "s.b"
        for (auto &sb : allBoards) {
            int s= sb.first, b= sb.second;
            double val = get32offset0(false,s,b);
            if (val<0) continue;

            int localSec= (s % 32);
            int phiMin= 8*localSec, phiMax= phiMin+7;
            int etaMin=0, etaMax=0;
            if (s<32) {
                etaMin=48 + 8*b;
                etaMax=etaMin+7;
            } else {
                int base= 40 - 8*b;
                etaMin= base; etaMax= base+7;
            }
            if (etaMin>etaMax) std::swap(etaMin,etaMax);

            double cx= 0.5*(phiMin + phiMax +1);
            double cy= 0.5*(etaMin + etaMax +1);

            TLatex lat;
            lat.SetTextAlign(22);
            lat.SetTextSize(0.015);
            lat.DrawLatex(cx,cy, Form("%d.%d",s,b));
        }

        // bottom => new(bestW)
        cH.cd(2);
        gPad->SetRightMargin(0.20);

        TH2F* h2Best = new TH2F("h2NewBest","NEW (Best W, #Delta=0);#phi index;#eta index",
                                256,0,256,96,0,96);
        h2Best->SetStats(false);
        h2Best->GetZaxis()->SetTitle("avg ADC");

        for (auto &sb : allBoards) {
            int s= sb.first, b= sb.second;
            double val= getBestOffset0(false, s, b);
            if(val<0) continue;

            int localSec= (s % 32);
            int phiMin= 8*localSec, phiMax= phiMin+7;
            int etaMin=0, etaMax=0;
            if (s<32) {
                etaMin=48 + 8*b;
                etaMax= etaMin+7;
            } else {
                int base=40 - 8*b;
                etaMin= base;
                etaMax= base+7;
            }
            if (etaMin>etaMax) std::swap(etaMin,etaMax);

            for(int xx=phiMin; xx<=phiMax; xx++){
                for(int yy=etaMin; yy<=etaMax; yy++){
                    h2Best->SetBinContent(xx+1, yy+1, val);
                }
            }
        }
        h2Best->Draw("COLZ");

        // label each IB
        for (auto &sb : allBoards) {
            int s= sb.first, b= sb.second;
            double val= getBestOffset0(false, s, b);
            if (val<0) continue;

            int localSec= (s % 32);
            int phiMin= 8*localSec, phiMax= phiMin+7;
            int etaMin=0, etaMax=0;
            if (s<32) {
                etaMin=48+8*b;
                etaMax=etaMin+7;
            } else {
                int base=40-8*b;
                etaMin= base;
                etaMax= base+7;
            }
            if (etaMin>etaMax) std::swap(etaMin,etaMax);

            double cx= 0.5*(phiMin + phiMax +1);
            double cy= 0.5*(etaMin + etaMax +1);

            TLatex lat;
            lat.SetTextAlign(22);
            lat.SetTextSize(0.015);
            lat.DrawLatex(cx,cy, Form("%d.%d",s,b));
        }

        std::string outC= "/Users/patsfan753/Desktop/PMTgainsAna/output/2Dmaps/heatMap_new_combined_colz.png";
        cH.SaveAs(outC.c_str());
        std::cout << "[INFO] => Created NEW sample combined heat maps: "<< outC <<"\n";

        delete h2W32;
        delete h2Best;
    }

    //----------------------------------------------------------------------
    // 9) The EXTRA STEP => 1D red + 2D highlight of top-3 spikes in newFitMap
    //----------------------------------------------------------------------
    {
        // (A) Build 1D distribution for “new bestW offset=0”
        TH1F *histNewBest = new TH1F("histNewBest","NEW bestW offset=0;ADC;Count",128,0,16384);
        histNewBest->SetLineColor(kRed);
        histNewBest->SetLineWidth(3);

        // binIndex => vector of (s,b)
        std::map<int, std::vector< std::pair<int,int> >> newBestContributors;

        // fill hist + record (s,b)
        for (auto &pr : newFitMap) {
            int s= pr.first.first;
            int b= pr.first.second;
            const BoardFitData &bfd = pr.second;

            int wBest = bfd.bestW;
            if (wBest < 1) continue; // invalid
            // find the vector of XY points for bestW
            auto itW = bfd.data.find(wBest);
            if (itW == bfd.data.end()) continue;

            for (auto &xy : itW->second) {
                if (std::fabs(xy.x)<1e-9 && xy.y>0.) {
                    double adcVal = xy.y;
                    int binIdx = histNewBest->FindBin(adcVal);
                    histNewBest->Fill(adcVal);
                    newBestContributors[binIdx].push_back({s,b});
                }
            }
        }

        // (B) Identify top-3 spikes in [2000..4000]
        int firstBin= histNewBest->FindBin(2000.);
        int lastBin = histNewBest->FindBin(4000.);
        std::vector< std::pair<int,int> > countsBins; // (count, binIdx)

        for(int ib= firstBin; ib<= lastBin; ib++){
            int c= (int)histNewBest->GetBinContent(ib);
            if(c>0) countsBins.push_back({c, ib});
        }
        // sort descending
        std::sort(countsBins.begin(), countsBins.end(),
                  [](auto &a, auto &b){ return a.first > b.first; });

        std::vector< std::tuple<int,int,double> > topBins; // (count, bin, xCenter)
        for(int i=0; i<3 && i<(int)countsBins.size(); i++){
            int count = countsBins[i].first;
            int bn    = countsBins[i].second;
            double xc = histNewBest->GetBinCenter(bn);
            topBins.push_back({count, bn, xc});
        }

        // gather all IBs from those top bins => highlight them
        std::vector<std::pair<int,int>> highlightIBs;
        for (auto &tb : topBins) {
            int binIdx = std::get<1>(tb);
            auto &vec  = newBestContributors[binIdx];
            highlightIBs.insert(highlightIBs.end(), vec.begin(), vec.end());
        }
        // remove duplicates
        std::sort(highlightIBs.begin(), highlightIBs.end());
        highlightIBs.erase(std::unique(highlightIBs.begin(),highlightIBs.end()),
                           highlightIBs.end());

        // (C) Build 2×1 canvas => top=1D red, bottom=2D color map
        TCanvas cX("cExtra","New BEST: 1D top, 2D bottom", 1200,1200);
        cX.Divide(1,2);

        // top => 1D
        cX.cd(1);
        gPad->SetRightMargin(0.05);

        double yMax= 1.2 * histNewBest->GetMaximum();
        histNewBest->SetMaximum(yMax);
        histNewBest->SetTitle("NEW sample: best TP (#Delta=0) => 1D distribution");
        histNewBest->GetXaxis()->SetTitle("Avg ADC");
        histNewBest->GetYaxis()->SetTitle("Count of IBs");
        histNewBest->Draw("HIST");

        // green lines => [2000..4000]
        TLine L1(2000,0,2000,yMax), L2(4000,0,4000,yMax);
        for (auto L : {&L1,&L2}) {
            L->SetLineColor(kGreen+2);
            L->SetLineStyle(2);
            L->SetLineWidth(2);
            L->Draw("same");
        }

        // black lines => top bin centers
        TLine peakDummy(0,0,0,0);
        peakDummy.SetLineColor(kBlack);
        peakDummy.SetLineStyle(2);
        peakDummy.SetLineWidth(2);

        for (auto &tb : topBins) {
            double xC= std::get<2>(tb);
            TLine *vLine= new TLine(xC, 0., xC, 0.9*yMax);
            vLine->SetLineColor(kBlack);
            vLine->SetLineStyle(2);
            vLine->SetLineWidth(2);
            vLine->Draw("same");
        }

        // small legend
        TLegend *leg = new TLegend(0.63, 0.62, 0.88, 0.75);
        leg->SetBorderSize(1);
        leg->SetFillColorAlpha(kWhite,0.6);
        leg->SetTextSize(0.03);
        leg->AddEntry(histNewBest, "NEW best TP(#Delta=0)","l");
        leg->AddEntry(&L1, "Target range: [2000..4000]","l");
        leg->AddEntry(&peakDummy,"Top-3 peak centers","l");
        leg->Draw();

        // text block => top bins
        {
            TLatex tx;
            tx.SetNDC(true);
            tx.SetTextAlign(13); // left
            tx.SetTextSize(0.03);

            double xT= 0.63, yT= 0.59, dy= 0.045;
            for (size_t i=0; i<topBins.size(); i++){
                int count = std::get<0>(topBins[i]);
                double xC = std::get<2>(topBins[i]);
                tx.DrawLatex(xT, yT - i*dy,
                             Form("ADC=%.0f => %d IBs", xC, count));
            }
        }

        // bottom => 2D color map
        cX.cd(2);
        gPad->SetRightMargin(0.20);

        TH2F* h2 = new TH2F("h2","NEW bestW offset=0 (2D);#phi;#eta",
                            256,0,256,96,0,96);
        h2->SetStats(false);
        h2->GetZaxis()->SetTitle("avg ADC");

        // fill from getBestOffset0(false => new)
        for (auto &pr : newFitMap) {
            int s= pr.first.first;
            int b= pr.first.second;
            double val= getBestOffset0(false, s,b);
            if (val<0) continue;

            int localSec= (s % 32);
            int phiMin= 8*localSec, phiMax= phiMin+7;
            int etaMin=0, etaMax=0;
            if(s<32) {
                etaMin=48+8*b;
                etaMax= etaMin+7;
            } else {
                etaMin=40-8*b;
                etaMax= etaMin+7;
            }
            if (etaMin>etaMax) std::swap(etaMin,etaMax);

            for(int xx= phiMin; xx<=phiMax; xx++){
                for(int yy= etaMin; yy<=etaMax; yy++){
                    h2->SetBinContent(xx+1, yy+1, val);
                }
            }
        }
        h2->Draw("COLZ");

        // label each IB => "s.b"
        for (auto &pr : newFitMap) {
            int s= pr.first.first;
            int b= pr.first.second;
            double val= getBestOffset0(false, s,b);
            if (val<0) continue;

            int localSec= (s % 32);
            int phiMin= 8*localSec, phiMax= phiMin+7;
            int etaMin=0, etaMax=0;
            if(s<32) {
                etaMin=48+8*b;
                etaMax= etaMin+7;
            } else {
                etaMin=40-8*b;
                etaMax= etaMin+7;
            }
            if(etaMin>etaMax) std::swap(etaMin,etaMax);

            double cx= 0.5*(phiMin + phiMax);
            double cy= 0.5*(etaMin + etaMax);

            TLatex lab;
            lab.SetTextAlign(22);
            lab.SetTextSize(0.015);
            lab.DrawLatex(cx,cy, Form("%d.%d", s,b));
        }

        // Over-plot open red circle for each highlight IB
        for (auto &sb : highlightIBs) {
            int s= sb.first;
            int b= sb.second;

            int localSec= (s % 32);
            double phiMin= 8*localSec, phiMax= phiMin+8;

            double etaMin, etaMax;
            if(s<32) {
                etaMin=48+8*b;
                etaMax= etaMin+8;
            } else {
                etaMin=40-8*b;
                etaMax= etaMin+8;
            }
            if (etaMin>etaMax) std::swap(etaMin,etaMax);

            double cx= 0.5*(phiMin + phiMax);
            double cy= 0.5*(etaMin + etaMax);

            TMarker *mk= new TMarker(cx, cy, 24);
            mk->SetMarkerColor(kRed);
            mk->SetMarkerSize(2.5);
            mk->Draw();
        }

        std::string outEx = "/Users/patsfan753/Desktop/PMTgainsAna/output/2Dmaps/newBest_1D_top_2D_bottom_highlight.png";
        cX.SaveAs(outEx.c_str());
        std::cout << "[INFO] => Created extra figure => " << outEx << "\n";

        delete h2;
        delete histNewBest;
    }

    //----------------------------------------------------------------------
    // 10) Print final tallies => old_total, new_total
    //----------------------------------------------------------------------
    int old_total = old_below + old_inrange + old_above + old_none;
    int new_total = new_below + new_inrange + new_above + new_none;

    std::cout << "\n[DEBUG] OLD(32,0) color-coded => total IB= "<< old_total << "\n"
              << "  below=  " << old_below <<" => "<< (100.*old_below/old_total)  <<"%\n"
              << "  inRange=" << old_inrange<<" => "<< (100.*old_inrange/old_total)<<"%\n"
              << "  above=  " << old_above  <<" => "<< (100.*old_above/old_total)  <<"%\n"
              << "  none=   " << old_none   <<" => "<< (100.*old_none/old_total)   <<"%\n";

    std::cout << "[DEBUG] NEW(32,0) color-coded => total IB= "<< new_total << "\n"
              << "  below=  " << new_below <<" => "<< (100.*new_below/new_total)  <<"%\n"
              << "  inRange=" << new_inrange<<" => "<< (100.*new_inrange/new_total)<<"%\n"
              << "  above=  " << new_above  <<" => "<< (100.*new_above/new_total)  <<"%\n"
              << "  none=   " << new_none   <<" => "<< (100.*new_none/new_total)   <<"%\n";
}



//--------------------------------------------------------------------------
//  RETURN:  std::pair<double,double>  ==> { mean_k_new , sigma_k_new }
//--------------------------------------------------------------------------
std::pair<double,double> produceKAndAmplitudePlots(
    const std::map<std::pair<int,int>, BoardFitData>& oldFitMap,
    const std::map<std::pair<int,int>, BoardFitData>& newFitMap)
{
    //----------------------------------------------------------------------
    // 0)  Output directory
    //----------------------------------------------------------------------
    const std::string outDir =
        "/Users/patsfan753/Desktop/PMTgainsAna/output/summary_ampAndk_vals";
    gSystem->Exec(Form("mkdir -p %s", outDir.c_str()));

    //----------------------------------------------------------------------
    // 1)  Collect parameters for boards present in BOTH samples
    //----------------------------------------------------------------------
    std::vector<Kval>   kPairs;     // (sector, ib, kOld, kNew)
    std::vector<AmpVal> aPairs;     // (sector, ib, ampOld, ampNew)

    for (const auto& pr : oldFitMap)
    {
        auto itN = newFitMap.find(pr.first);
        if (itN == newFitMap.end()) continue;  // board not in new => skip

        kPairs.push_back({
            pr.first.first,
            pr.first.second,
            pr.second.kParam,
            itN->second.kParam
        });

        aPairs.push_back({
            pr.first.first,
            pr.first.second,
            pr.second.amplitude,
            itN->second.amplitude
        });
    }

    if (kPairs.empty())
    {
        std::cout << "[WARN] → No common boards between OLD & NEW – skipping plots.\n";
        return {0.,0.};
    }

    //----------------------------------------------------------------------
    // 2)  k-value histograms (OLD vs NEW) => PNG #1
    //----------------------------------------------------------------------
    double kMax = 0.;
    for (auto &kp : kPairs) {
        kMax = std::max(kMax, std::max(kp.kOld, kp.kNew));
    }
    kMax = (kMax < 1e-7) ? 1e-5 : 1.2 * kMax;

    // Build histograms for k
    TH1F hKold("hKold","Old k; k; # IBs", 50,0.,kMax);
    TH1F hKnew("hKnew","New k; k; # IBs", 50,0.,kMax);
    hKold.SetLineColor(kBlue+2);
    hKold.SetLineWidth(3);
    hKnew.SetLineColor(kRed+1);
    hKnew.SetLineWidth(3);

    // Fill histograms
    for (auto &kp : kPairs) {
        if (kp.kOld>0) hKold.Fill(kp.kOld);
        if (kp.kNew>0) hKnew.Fill(kp.kNew);
    }

    TCanvas cK("cK","k distributions",1200,600);
    cK.Divide(2,1);

    // --- LEFT side => "Old k" ---
    cK.cd(1);
    hKold.Draw("HIST");

    // Step 1: quick fit across entire range
    TF1 fG_old_init("fG_old_init","gaus",0.,kMax);
    fG_old_init.SetParameters(hKold.GetMaximum(), hKold.GetMean(), hKold.GetRMS());
    TFitResultPtr frOld1 = hKold.Fit(&fG_old_init,"Q S R");
    
    // from that, define narrower range
    double meanO1     = frOld1->Parameter(1);
    double sigmaO1    = std::fabs(frOld1->Parameter(2));
    double loOld      = std::max(0.0, meanO1 - 3.*sigmaO1);
    double hiOld      = meanO1 + 3.*sigmaO1;

    // Step 2: final fit in narrower range
    TF1 fG_old("fG_old","gaus",loOld,hiOld);
    fG_old.SetLineColor(kGreen+3);
    fG_old.SetLineWidth(2);

    // set initial guesses from previous step
    fG_old.SetParameters(frOld1->Parameter(0), meanO1, sigmaO1);
    TFitResultPtr frOld2 = hKold.Fit(&fG_old,"Q S R+");

    // Retrieve final fit parameters + errors
    double mu_old       = frOld2->Parameter(1);
    double sigma_old    = frOld2->Parameter(2);
    double emu_old      = frOld2->ParError(1);
    double esigma_old   = frOld2->ParError(2);
    double chi2_old     = frOld2->Chi2();
    double ndf_old      = frOld2->Ndf();
    double chi2NdfOld   = (ndf_old>0 ? chi2_old/ndf_old : 0.0);
    double resolutionOld= (std::fabs(mu_old)>1e-15 ? 100.*sigma_old/std::fabs(mu_old) : 0.);

    fG_old.Draw("SAME");

    // Print fit info with TLatex
    {
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextSize(0.04);

        double y=0.85;
        lat.DrawLatex(0.55,y, Form("#mu = %.4g #pm %.4g", mu_old, emu_old));
        y -=0.06;
        lat.DrawLatex(0.55,y, Form("#sigma = %.4g #pm %.4g", sigma_old, esigma_old));
        y -=0.06;
        lat.DrawLatex(0.55,y, Form("#chi^{2}/NDF = %.2f", chi2NdfOld));
        y -=0.06;
        lat.DrawLatex(0.55,y, Form("Res=%.1f%%", resolutionOld));
        y -=0.06;
        lat.DrawLatex(0.55,y, Form("Entries=%.0f", hKold.GetEntries()));
    }

    // --- RIGHT side => "New k" ---
    cK.cd(2);
    hKnew.Draw("HIST");

    // Step 1: initial fit on entire range
    TF1 fG_new_init("fG_new_init","gaus",0.,kMax);
    fG_new_init.SetParameters(hKnew.GetMaximum(), hKnew.GetMean(), hKnew.GetRMS());
    TFitResultPtr frNew1 = hKnew.Fit(&fG_new_init,"Q S R");

    double meanN1     = frNew1->Parameter(1);
    double sigmaN1    = std::fabs(frNew1->Parameter(2));
    double loNew      = std::max(0.0, meanN1 - 3.*sigmaN1);
    double hiNew      = meanN1 + 3.*sigmaN1;

    // Step 2: final narrower fit
    TF1 fG_new("fG_new","gaus",loNew,hiNew);
    fG_new.SetLineColor(kGreen+3);
    fG_new.SetLineWidth(2);
    fG_new.SetParameters(frNew1->Parameter(0), meanN1, sigmaN1);

    TFitResultPtr frNew2 = hKnew.Fit(&fG_new,"Q S R+");

    double mu_new       = frNew2->Parameter(1);
    double sigma_new    = frNew2->Parameter(2);
    double emu_new      = frNew2->ParError(1);
    double esigma_new   = frNew2->ParError(2);
    double chi2_new     = frNew2->Chi2();
    double ndf_new      = frNew2->Ndf();
    double chi2NdfNew   = (ndf_new>0 ? chi2_new/ndf_new : 0.);
    double resolutionNew= (std::fabs(mu_new)>1e-15 ? 100.*sigma_new/std::fabs(mu_new) : 0.);

    fG_new.Draw("SAME");

    // Print fit info with TLatex
    {
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextSize(0.04);

        double y=0.85;
        lat.DrawLatex(0.55,y, Form("#mu = %.4g #pm %.4g", mu_new, emu_new));
        y -=0.06;
        lat.DrawLatex(0.55,y, Form("#sigma = %.4g #pm %.4g", sigma_new, esigma_new));
        y -=0.06;
        lat.DrawLatex(0.55,y, Form("#chi^{2}/NDF = %.2f", chi2NdfNew));
        y -=0.06;
        lat.DrawLatex(0.55,y, Form("Res=%.1f%%", resolutionNew));
        y -=0.06;
        lat.DrawLatex(0.55,y, Form("Entries=%.0f", hKnew.GetEntries()));
    }

    // We'll return these final new k-values
    double muK_new    = mu_new;
    double sigmaK_new = sigma_new;

    cK.SaveAs( (outDir + "/kValue_distributions.png").c_str() );
    std::cout << "[INFO] → wrote kValue_distributions.png\n";

    //----------------------------------------------------------------------
    // 3)  amplitude histograms (Gaussian) => PNG #2
    //----------------------------------------------------------------------
    // Use separate X ranges for old vs new to avoid wasted space.
    double aMaxOld = 0.;
    double aMaxNew = 0.;
    for (auto &a : aPairs) {
        if (a.ampOld>0 && a.ampOld>aMaxOld) aMaxOld = a.ampOld;
        if (a.ampNew>0 && a.ampNew>aMaxNew) aMaxNew = a.ampNew;
    }
    if (aMaxOld<1e-7) aMaxOld= 1e-5; else aMaxOld *= 1.2;
    if (aMaxNew<1e-7) aMaxNew= 1e-5; else aMaxNew *= 1.2;

    TH1F hAold("hAold","Old amplitude;Amplitude;# IBs",50, 0., aMaxOld);
    TH1F hAnew("hAnew","New amplitude;Amplitude;# IBs",50, 0., aMaxNew);

    hAold.SetLineColor(kBlue+2);
    hAold.SetLineWidth(3);
    hAnew.SetLineColor(kRed+1);
    hAnew.SetLineWidth(3);

    // Fill amplitude hists
    for (auto &a : aPairs) {
        if (a.ampOld>0) hAold.Fill(a.ampOld);
        if (a.ampNew>0) hAnew.Fill(a.ampNew);
    }

    TCanvas cA("cA","Amplitude distributions",1200,600);
    cA.Divide(2,1);

    // =============== Left => old amplitude (Gaussian fit) ===============
    cA.cd(1);
    hAold.Draw("HIST");

    TF1 fG_aold_init("fG_aold_init","gaus", 0., aMaxOld);
    fG_aold_init.SetParameters(hAold.GetMaximum(), hAold.GetMean(), hAold.GetRMS());
    TFitResultPtr frOldAmp1 = hAold.Fit(&fG_aold_init,"Q S R");

    double muO1Amp     = frOldAmp1->Parameter(1);
    double sigmaO1Amp  = std::fabs(frOldAmp1->Parameter(2));
    double loOldAmp    = std::max(0.0, muO1Amp - 3.*sigmaO1Amp);
    double hiOldAmp    = std::min(aMaxOld, muO1Amp + 3.*sigmaO1Amp);

    TF1 fG_aold("fG_aold","gaus", loOldAmp, hiOldAmp);
    fG_aold.SetLineColor(kGreen+3);
    fG_aold.SetLineWidth(2);
    fG_aold.SetParameters(frOldAmp1->Parameter(0), muO1Amp, sigmaO1Amp);

    TFitResultPtr frOldAmp2 = hAold.Fit(&fG_aold,"Q S R+");

    double muOldAmp     = frOldAmp2->Parameter(1);
    double sigmaOldAmp  = frOldAmp2->Parameter(2);
    double emuOldAmp    = frOldAmp2->ParError(1);
    double esigOldAmp   = frOldAmp2->ParError(2);

    double chi2OldAmp   = frOldAmp2->Chi2();
    double ndfOldAmp    = frOldAmp2->Ndf();
    double chi2NdfOldAmp= (ndfOldAmp>0 ? chi2OldAmp/ndfOldAmp : 0.);
    double resOldAmp    = (std::fabs(muOldAmp)>1e-12 ? 100.*sigmaOldAmp/std::fabs(muOldAmp) : 0.);

    fG_aold.Draw("SAME");

    // TLatex annotation
    {
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextSize(0.04);
        double y=0.85;

        lat.DrawLatex(0.55,y, Form("#mu = %.1f #pm %.1f", muOldAmp, emuOldAmp)); y-=0.06;
        lat.DrawLatex(0.55,y, Form("#sigma = %.1f #pm %.1f", sigmaOldAmp, esigOldAmp)); y-=0.06;
        lat.DrawLatex(0.55,y, Form("#chi^{2}/NDF = %.2f", chi2NdfOldAmp));         y-=0.06;
        lat.DrawLatex(0.55,y, Form("Res=%.1f%%", resOldAmp));                     y-=0.06;
        lat.DrawLatex(0.55,y, Form("Entries=%.0f", hAold.GetEntries()));
    }

    // =============== Right => new amplitude (Gaussian fit) ===============
    cA.cd(2);
    hAnew.Draw("HIST");

    TF1 fG_anew_init("fG_anew_init","gaus",0., aMaxNew);
    fG_anew_init.SetParameters(hAnew.GetMaximum(), hAnew.GetMean(), hAnew.GetRMS());
    TFitResultPtr frNewAmp1 = hAnew.Fit(&fG_anew_init,"Q S R");

    double muN1Amp     = frNewAmp1->Parameter(1);
    double sigmaN1Amp  = std::fabs(frNewAmp1->Parameter(2));
    double loNewAmp    = std::max(0.0, muN1Amp - 3.*sigmaN1Amp);
    double hiNewAmp    = std::min(aMaxNew, muN1Amp + 3.*sigmaN1Amp);

    TF1 fG_anew("fG_anew","gaus",loNewAmp,hiNewAmp);
    fG_anew.SetLineColor(kGreen+3);
    fG_anew.SetLineWidth(2);
    fG_anew.SetParameters(frNewAmp1->Parameter(0), muN1Amp, sigmaN1Amp);

    TFitResultPtr frNewAmp2 = hAnew.Fit(&fG_anew,"Q S R+");

    double muNewAmp    = frNewAmp2->Parameter(1);
    double sigmaNewAmp = frNewAmp2->Parameter(2);
    double emuNewAmp   = frNewAmp2->ParError(1);
    double esigNewAmp  = frNewAmp2->ParError(2);

    double chi2NewAmp   = frNewAmp2->Chi2();
    double ndfNewAmp    = frNewAmp2->Ndf();
    double chi2NdfNewAmp= (ndfNewAmp>0 ? chi2NewAmp/ndfNewAmp : 0.);
    double resNewAmp    = (std::fabs(muNewAmp)>1e-12 ? 100.*sigmaNewAmp/std::fabs(muNewAmp) : 0.);

    fG_anew.Draw("SAME");

    // TLatex annotation
    {
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextSize(0.04);
        double y=0.85;

        lat.DrawLatex(0.55,y, Form("#mu = %.1f #pm %.1f", muNewAmp, emuNewAmp)); y-=0.06;
        lat.DrawLatex(0.55,y, Form("#sigma = %.1f #pm %.1f", sigmaNewAmp, esigNewAmp)); y-=0.06;
        lat.DrawLatex(0.55,y, Form("#chi^{2}/NDF = %.2f", chi2NdfNewAmp));            y-=0.06;
        lat.DrawLatex(0.55,y, Form("Res=%.1f%%", resNewAmp));                         y-=0.06;
        lat.DrawLatex(0.55,y, Form("Entries=%.0f", hAnew.GetEntries()));
    }

    // Save => amplitude_distributions_gaus.png
    cA.SaveAs( (outDir + "/amplitude_distributions_gaus.png").c_str() );
    std::cout << "[INFO] → wrote amplitude_distributions_gaus.png\n";

    //----------------------------------------------------------------------
    // 4)  k-ratio NEW/OLD => NO FIT, show points => PNG #3
    //----------------------------------------------------------------------
    std::vector<double> kRatio;
    double kRmax = 0.;

    // We'll track which boards exceed ratio=10 for printing
    std::vector< std::pair<std::pair<int,int>, double> > outliers;

    for (auto &kp : kPairs)
    {
        if (kp.kOld>0 && kp.kNew>0) {
            double r = kp.kNew / kp.kOld;
            if(r>10.0) {
                outliers.push_back({{kp.sector, kp.ib}, r});
            }
            else {
                kRatio.push_back(r);
            }
            if (r>kRmax) kRmax = r;
        }
    }

    if (!kRatio.empty())
    {
        // fix range => 0..10
        TH1F hKR("hKR","k ratio (new/old);Ratio;# IBs",50, 0., 10.);
        for (double r : kRatio) {
            hKR.Fill(r);
        }

        TCanvas cKR("cKR","k ratio",900,600);

        hKR.SetMarkerStyle(20);
        hKR.SetMarkerColor(kBlue+2);
        hKR.SetLineColor(kBlue+2);
        hKR.SetLineWidth(2);

        hKR.Draw("PE");

        // Mark outliers
        if(!outliers.empty()) {
            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextSize(0.03);

            double yPos = 0.8;
            lat.DrawLatex(0.15,yPos,"#bf{IBs OUT OF RANGE (ratio>10):}");
            yPos -= 0.05;

            int countPrinted=0;
            for(auto &ol : outliers) {
                if(countPrinted>=6) break;
                lat.DrawLatex(
                    0.15,yPos,
                    Form("S=%d IB=%d => ratio=%.2f",
                         ol.first.first,
                         ol.first.second,
                         ol.second));
                yPos -= 0.05;
                countPrinted++;
            }
        }

        cKR.SaveAs( (outDir + "/kRatio_bestW_gaus.png").c_str() );
        std::cout << "[INFO] → wrote kRatio_bestW_gaus.png\n";
    }
    else {
        std::cout << "[WARN] → kRatio histogram skipped (no valid data).\n";
    }

    //----------------------------------------------------------------------
    // 5)  amplitude-ratio NEW/OLD => NO FIT, show points => PNG #4
    //----------------------------------------------------------------------
    std::vector<double> aRatio;
    double aRmax = 0.;
    for (auto &av : aPairs)
    {
        if (av.ampOld>0 && av.ampNew>0) {
            double r = av.ampNew/av.ampOld;
            aRatio.push_back(r);
            if (r>aRmax) aRmax=r;
        }
    }

    if (!aRatio.empty())
    {
        aRmax *= 1.2;
        TH1F hAR("hAR","Amplitude ratio (new/old);Ratio;# IBs",50,0.,aRmax);
        for (double r : aRatio) {
            hAR.Fill(r);
        }

        TCanvas cAR("cAR","Amplitude ratio",900,600);

        hAR.SetMarkerStyle(20);
        hAR.SetMarkerColor(kBlue+2);
        hAR.SetLineColor(kBlue+2);
        hAR.SetLineWidth(2);

        hAR.Draw("PE");

        cAR.SaveAs( (outDir + "/amplitudeRatio_bestW_gaus.png").c_str() );
        std::cout << "[INFO] → wrote amplitudeRatio_bestW_gaus.png\n";
    }
    else {
        std::cout << "[WARN] → amplitudeRatio histogram skipped (no valid data).\n";
    }

    //----------------------------------------------------------------------
    // 6) Done
    //----------------------------------------------------------------------
    std::cout << "[INFO] → produceKAndAmplitudePlots() finished (4 PNGs written).\n\n";

    // Return final new k fit parameters from step #2
    return { muK_new , sigmaK_new };
}
// -------------------------------------------------------------------------
//  Classify NEW boards whose k-parameter is > 1 σ away from the global
//  NEW-sample ⟨k⟩ and compare their BEST-W amplitudes.
//
//  Parameters
//  ----------
//   newFitMap   : results for NEW sample (key = {sector,ib})
//   oldFitMap   : results for OLD sample   (only needed for the PW summary)
//   meanK_new   : ⟨k⟩  from the NEW-sample Gaussian fit
//   sigmaK_new  : σ(k) from the NEW-sample Gaussian fit
//   outputPNG   : full pathname for the 2-panel PNG that will be produced
// -------------------------------------------------------------------------
// -----------------------------------------------------------------------------
//  Classify NEW boards whose k-parameter is > 1 σ away from the global
//  NEW-sample ⟨k⟩ and compare their BEST-W amplitudes.
//
//  Then for the “in-sigma” amplitude distribution (|k - meanK_new| < sigmaK_new)
//  and the “out-sigma” amplitude distribution (|k - meanK_new| > sigmaK_new),
//  we fit each amplitude histogram with a Gaussian, and print the IBs that
//  are ±3σ away from that fitted mean amplitude (showing each board's amplitude
//  and kParam on the top-left corner).
// -----------------------------------------------------------------------------
void outOfSigmaAmplitudeCheck(
        const std::map< std::pair<int,int>, BoardFitData > &newFitMap,
        const std::map< std::pair<int,int>, BoardFitData > &oldFitMap,
        double  meanK_new,
        double  sigmaK_new,
        const   std::string &outputPNG)
{
    if (sigmaK_new < 1e-12) sigmaK_new = 1e-12;  // guard against near-zero σ

    std::vector<double> amps_inSigma, amps_outSigma;

    // We also keep track of the actual boards so we can do the 3σ amplitude check
    // after we do the Gaussian fit:
    //   - inSigmaBoards => all boards { (s,ib), amplitude, k }
    //   - outSigmaBoards => likewise
    struct BoardAmpK {
        std::pair<int,int> sb; // (sector, ib)
        double amp;
        double k;
    };
    std::vector<BoardAmpK> inSigmaBoards, outSigmaBoards;

    //----------------------------------------------------------------------
    // 1)  Scan NEW boards once — decide “in / out” and collect amplitudes
    //----------------------------------------------------------------------
    for (const auto &pr : newFitMap)
    {
        int s = pr.first.first;
        int b = pr.first.second;

        const BoardFitData &bfd = pr.second;
        if (bfd.kParam <= 0.) continue;         // skip invalid k

        double amp = bfd.amplitude;             // amplitude from bestW
        if (amp <= 0.) continue;                // skip invalid A

        // This is the “k-based” cut from original code
        bool isOut = (std::fabs(bfd.kParam - meanK_new) > sigmaK_new);

        if (isOut) {
            amps_outSigma.push_back(amp);
            outSigmaBoards.push_back({{s,b}, amp, bfd.kParam});

            std::cout << "  [OUTLIER by k-sigma]  s=" << s
                      << "  ib="       << b
                      << "  k="        << bfd.kParam
                      << "  A="        << amp << '\n';
        }
        else {
            amps_inSigma.push_back(amp);
            inSigmaBoards.push_back({{s,b}, amp, bfd.kParam});
        }
    }

    std::cout << "\n[CHECK]  # in-σ boards  = " << amps_inSigma.size()
              << "\n         # out-σ boards = " << amps_outSigma.size() << "\n";

    //----------------------------------------------------------------------
    // 2)  Build amplitude histograms
    //----------------------------------------------------------------------
    double maxA = 0.;
    for (double v : amps_inSigma ) if (v > maxA) maxA = v;
    for (double v : amps_outSigma) if (v > maxA) maxA = v;
    if (maxA < 1e-7) maxA = 1e-5;

    TH1F hA_in ("hA_in" ,
                "Amplitude (|k-#LTk#GT|<#sigma);Amplitude;Count",
                50, 0., 1.2*maxA);
    TH1F hA_out("hA_out",
                "Amplitude (|k-#LTk#GT|>#sigma);Amplitude;Count",
                50, 0., 1.2*maxA);

    for (double v : amps_inSigma ) hA_in .Fill(v);
    for (double v : amps_outSigma) hA_out.Fill(v);

    //----------------------------------------------------------------------
    // 3)  Simple stats on amplitude arrays
    //----------------------------------------------------------------------
    auto meanVec  = [](const std::vector<double>& v)
                    { return v.empty()? 0. : TMath::Mean (v.size(), &v[0]); };
    auto medianVec= [](const std::vector<double>& v)
                    { return v.empty()? 0. : TMath::Median(v.size(), &v[0]); };

    double meanIn    = meanVec(amps_inSigma);
    double medianIn  = medianVec(amps_inSigma);
    double meanOut   = meanVec(amps_outSigma);
    double medianOut = medianVec(amps_outSigma);

    //----------------------------------------------------------------------
    // 4) Two-panel canvas => we add Gaussian fits on each amplitude hist
    //    then find boards that are ±3σ from each fitted mean
    //----------------------------------------------------------------------
    TCanvas cChk("cChk","Amplitude check vs k-outliers",1200,600);
    cChk.Divide(2,1);

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.03);

    // We define function objects for each hist => gaus
    TF1 fG_in("fG_in", "gaus", 0., 1.2*maxA);
    TF1 fG_out("fG_out","gaus", 0., 1.2*maxA);

    // The sets of boards (sector, ib, amplitude, k) that are ±3σ away from
    // the amplitude distribution's fitted mean
    std::vector<BoardAmpK> outliersIn, outliersOut;

    //-------------------------- in-σ panel --------------------------------
    cChk.cd(1);
    hA_in.SetLineColor(kBlue);
    hA_in.SetLineWidth(2);
    hA_in.Draw("HIST");

    if (hA_in.GetEntries() > 10) {
        hA_in.Fit(&fG_in, "Q R");
        double mu_in    = fG_in.GetParameter(1);
        double sigma_in = fG_in.GetParameter(2);

        // Mark boards that are ±3σ from mu_in
        for (auto &board : inSigmaBoards) {
            double diff = std::fabs(board.amp - mu_in);
            if (diff > 3.*sigma_in) {
                outliersIn.push_back(board);
            }
        }

        // Print the fitted result
        lat.DrawLatex(0.57,0.80,Form("#LTk#GT=%.3g, #sigma(k)=%.3g",
                                     meanK_new, sigmaK_new));
        lat.DrawLatex(0.57,0.74,Form("Fit #mu_{A}=%.1f, #sigma_{A}=%.1f", mu_in,sigma_in));
    }
    else {
        lat.DrawLatex(0.57,0.80,"(Not enough stats to fit in-sigma hist)");
    }

    lat.DrawLatex(0.57,0.66,Form("N_{IB} = %zu", amps_inSigma.size()));
    lat.DrawLatex(0.57,0.58,Form("#LT A #GT = %.1f", meanIn));
    lat.DrawLatex(0.57,0.52,Form("median(A)= %.1f", medianIn));

    // We'll print the outlier boards near the left side (just a few).
    {
        double yText = 0.40;
        int countShow=0;
        if(!outliersIn.empty()){
            lat.DrawLatex(0.15,yText,"#bf{±3#sigma Boards:}");
            yText -= 0.05;
        }
        for (auto &bo : outliersIn) {
            if(countShow>=5) break;
            lat.DrawLatex(0.15,yText,
                          Form("S=%d IB=%d A=%.1f k=%.3g",
                               bo.sb.first, bo.sb.second, bo.amp, bo.k));
            yText -= 0.05;
            countShow++;
        }
    }

    //------------------------- out-σ panel ---------------------------------
    cChk.cd(2);
    hA_out.SetLineColor(kRed);
    hA_out.SetLineWidth(2);
    hA_out.Draw("HIST");

    if (hA_out.GetEntries() > 10) {
        hA_out.Fit(&fG_out, "Q R");
        double mu_out    = fG_out.GetParameter(1);
        double sigma_out = fG_out.GetParameter(2);

        // Mark boards that are ±3σ from mu_out
        for (auto &board : outSigmaBoards) {
            double diff = std::fabs(board.amp - mu_out);
            if (diff > 3.*sigma_out) {
                outliersOut.push_back(board);
            }
        }

        // Print the fitted result
        lat.DrawLatex(0.57,0.80,Form("#LTk#GT=%.3g, #sigma(k)=%.3g",
                                     meanK_new, sigmaK_new));
        lat.DrawLatex(0.57,0.74,Form("Fit #mu_{A}=%.1f, #sigma_{A}=%.1f", mu_out,sigma_out));
    }
    else {
        lat.DrawLatex(0.57,0.80,"(Not enough stats to fit out-sigma hist)");
    }

    lat.DrawLatex(0.57,0.66,Form("N_{IB} = %zu", amps_outSigma.size()));
    lat.DrawLatex(0.57,0.58,Form("#LT A #GT = %.1f", meanOut));
    lat.DrawLatex(0.57,0.52,Form("median(A)= %.1f", medianOut));

    // We'll print the outlier boards near the left side (just a few).
    {
        double yText = 0.40;
        int countShow=0;
        if(!outliersOut.empty()){
            lat.DrawLatex(0.15,yText,"#bf{±3#sigma Boards:}");
            yText -= 0.05;
        }
        for (auto &bo : outliersOut) {
            if(countShow>=5) break;
            lat.DrawLatex(0.15,yText,
                          Form("S=%d IB=%d A=%.1f k=%.3g",
                               bo.sb.first, bo.sb.second, bo.amp, bo.k));
            yText -= 0.05;
            countShow++;
        }
    }

    cChk.SaveAs(outputPNG.c_str());
    std::cout << "[CHECK] → wrote " << outputPNG << "\n";

    //----------------------------------------------------------------------
    // 5)  (unchanged)  Best-pulse-width usage table 24–40 ns
    //----------------------------------------------------------------------
    {
        // We'll keep this exactly as in your original code
        std::map<int,int> cntOld, cntNew;
        for (int w=24; w<=40; ++w) { cntOld[w]=0; cntNew[w]=0; }

        for (auto &pr : oldFitMap) {
            int bw = pr.second.bestW;
            if (bw>=24 && bw<=40) {
                ++cntOld[bw];
            }
        }
        for (auto &pr : newFitMap) {
            int bw = pr.second.bestW;
            if (bw>=24 && bw<=40) {
                ++cntNew[bw];
            }
        }

        const char* G="\033[32m";
        const char* R="\033[31m";
        std::cout << "\033[93m=== BEST-W usage 24-40 ns ===\033[0m\n";
        std::cout << " Width |  OLD |  NEW\n";
        for (int w=24; w<=40; ++w)
            std::cout << "  " << std::setw(3) << w << "  |  "
                      << ((cntOld[w]?G:R)) << std::setw(3) << cntOld[w] << "\033[0m |  "
                      << ((cntNew[w]?G:R)) << std::setw(3) << cntNew[w] << "\033[0m\n";
        std::cout << "================================\n\n";
    }

    std::cout << "[CHECK] → out-of-σ amplitude analysis done.\n";
}


/************************  fixed, deterministic colour palette  *************/
static const std::vector<Color_t> kPalette = {
    kBlue  , kRed  , kGreen+2, kMagenta,
    kCyan+2, kOrange+7, kViolet,  kGray+2,
    kSpring+5, kAzure+7, kPink+7, kTeal+2,
    kOrange-3, kYellow+2, kBlue-6, kRed-6   // add more if you have >16 widths
};

/**
 * We keep a single static `knownWs` that is filled with *all* widths
 * found in oldFitMap and newFitMap the first time `produceMultiWOverlayPlots` is called.
 * Afterwards, `colorForWidth(w)` will always use that same ordering.
 */
static std::vector<int> gAllKnownWs;  // empty initially

// A helper that returns a color index for a given pulse width, based on gAllKnownWs.
Color_t colorForWidth(int w)
{
    // Insert w into gAllKnownWs if not already present
    // This ensures that if we ever see a new width, we keep it in the global list.
    if(std::find(gAllKnownWs.begin(), gAllKnownWs.end(), w) == gAllKnownWs.end()){
        gAllKnownWs.push_back(w);
        std::sort(gAllKnownWs.begin(), gAllKnownWs.end());
    }

    // Find index of w in gAllKnownWs
    auto it = std::find(gAllKnownWs.begin(), gAllKnownWs.end(), w);
    int idx = (int)std::distance(gAllKnownWs.begin(), it);

    // Wrap around if > palette size
    int cIndex = idx % (int)kPalette.size();
    return kPalette[cIndex];
}

// A small helper to get the overall (minY, maxY) from a TMultiGraph’s data
std::pair<double,double> getMinMaxY(const TMultiGraph* mg)
{
    double minY=1e15, maxY=-1e15;
    if(!mg || !mg->GetListOfGraphs()) return {0., 1.};

    for(int i=0; i< mg->GetListOfGraphs()->GetSize(); i++){
        auto obj = mg->GetListOfGraphs()->At(i);
        auto grE = dynamic_cast<TGraphErrors*>(obj);
        if(!grE) continue;
        for(int p=0; p<grE->GetN(); p++){
            double x=0,y=0;
            grE->GetPoint(p,x,y);
            if(y < minY) minY=y;
            if(y > maxY) maxY=y;
        }
    }
    if(minY>maxY) { minY=0.; maxY=1.; } // guard
    return {minY, maxY};
}

void produceMultiWOverlayPlots(
    const std::map< std::pair<int,int>, BoardFitData > &oldFitMap,
    const std::map< std::pair<int,int>, BoardFitData > &newFitMap,
    const std::string &multiW,       // e.g. outDir + "/multiW"
    const std::string &summary,      // e.g. outDir + "/summary"
    std::vector<Kval> &kRatios,      // OUTPUT => fill these
    std::vector<AmpVal> &amplitudeRatios // OUTPUT => fill these
)
{
    std::cout << "[DEBUG] Step5 => multi-W overlay plots...\n";

    // ------------------------------------------------------------------------
    //  STEP A: gather (union) of all widths from old and new data to fix color
    // ------------------------------------------------------------------------
    static bool firstCall = true;
    if (firstCall) {
        // Add known widths from oldFitMap
        for(const auto &pr : oldFitMap){
            const auto &bfd = pr.second;
            for(const auto &wk : bfd.data){
                int w = wk.first;
                if(std::find(gAllKnownWs.begin(), gAllKnownWs.end(), w) == gAllKnownWs.end())
                    gAllKnownWs.push_back(w);
            }
        }
        // Add known widths from newFitMap
        for(const auto &pr : newFitMap){
            const auto &bfd = pr.second;
            for(const auto &wk : bfd.data){
                int w = wk.first;
                if(std::find(gAllKnownWs.begin(), gAllKnownWs.end(), w) == gAllKnownWs.end())
                    gAllKnownWs.push_back(w);
            }
        }
        // Sort them once
        std::sort(gAllKnownWs.begin(), gAllKnownWs.end());
        firstCall = false;
    }
    // By now, gAllKnownWs is the union of old & new widths, in ascending order.

    int validIB=0;
    // We pick which widths to overlay. For example, your original {26..32}:
    // If you do want more/less widths, just adjust the wSet here.
    std::set<int> wSet{26,27,28,29,30,31,32};
    
    // For each board in oldFitMap that also exists in newFitMap:
    for (auto &oldPair : oldFitMap)
    {
        int s = oldPair.first.first;
        int b = oldPair.first.second;
        
        auto itNew = newFitMap.find({s,b});
        if (itNew == newFitMap.end()) continue;  // skip if not in new

        const BoardFitData &bfdOld = oldPair.second;
        const BoardFitData &bfdNew = itNew->second;
        if (bfdOld.data.empty() && bfdNew.data.empty()) continue;
        
        validIB++;
        
        // TMultiGraphs: mgOld, mgNew, mgRatio
        TMultiGraph* mgOld   = new TMultiGraph();
        mgOld->SetTitle(Form("OLD (S%d,IB%d);#Delta(mV);AvgADC", s,b));
        
        TMultiGraph* mgNew   = new TMultiGraph();
        mgNew->SetTitle(Form("NEW (S%d,IB%d);#Delta(mV);AvgADC", s,b));
        
        TMultiGraph* mgRatio = new TMultiGraph();
        mgRatio->SetTitle("New/Old;#Delta(mV);Ratio");

        // Add TGraphErrors for each width
        for (int w : wSet)
        {
            auto itOldW = bfdOld.data.find(w);
            auto itNewW = bfdNew.data.find(w);

            const std::vector<XYerr>* vecO =
                (itOldW!=bfdOld.data.end() ? &itOldW->second : nullptr);
            const std::vector<XYerr>* vecN =
                (itNewW!=bfdNew.data.end() ? &itNewW->second : nullptr);

            // color from the (now global) palette
            Color_t thisColor = colorForWidth(w);

            // OLD => marker style 20
            if(vecO && !vecO->empty()){
                TGraphErrors* grO = new TGraphErrors((int)vecO->size());
                for(int i=0;i<(int)vecO->size();i++){
                    grO->SetPoint(i, (*vecO)[i].x, (*vecO)[i].y);
                    grO->SetPointError(i, 0., (*vecO)[i].e);
                }
                grO->SetMarkerStyle(20);
                grO->SetMarkerColor(thisColor);
                grO->SetLineColor(thisColor);
                grO->SetNameTitle(Form("oldW%d",w), Form("Old(W=%d)",w));
                mgOld->Add(grO, "P");
            }

            // NEW => marker style 21
            if(vecN && !vecN->empty()){
                TGraphErrors* grN = new TGraphErrors((int)vecN->size());
                for(int i=0;i<(int)vecN->size();i++){
                    grN->SetPoint(i, (*vecN)[i].x, (*vecN)[i].y);
                    grN->SetPointError(i, 0., (*vecN)[i].e);
                }
                grN->SetMarkerStyle(21);
                grN->SetMarkerColor(thisColor);
                grN->SetLineColor(thisColor);
                grN->SetNameTitle(Form("newW%d",w), Form("New(W=%d)",w));
                mgNew->Add(grN, "P");
            }

            // RATIO => marker style 25
            if(vecO && vecN && !vecO->empty() && !vecN->empty()){
                // build a map of old x->y so we can compute ratio
                std::map<double,double> oldMap;
                for (auto &ptO : *vecO) {
                    oldMap[ptO.x] = ptO.y;
                }
                // ratio TGraph
                std::vector<double> rx, ry;
                for (auto &ptN : *vecN){
                    double dx = ptN.x;
                    if(oldMap.count(dx)>0 && oldMap[dx]>0){
                        rx.push_back(dx);
                        ry.push_back(ptN.y / oldMap[dx]);
                    }
                }
                if(!rx.empty()){
                    TGraph* grR = new TGraph((int)rx.size());
                    for(int i=0;i<(int)rx.size();i++){
                        grR->SetPoint(i, rx[i], ry[i]);
                    }
                    grR->SetMarkerStyle(25);
                    grR->SetMarkerColor(thisColor);
                    grR->SetLineColor(thisColor);
                    grR->SetNameTitle(Form("ratioW%d", w),
                                      Form("Ratio(W=%d)", w));
                    mgRatio->Add(grR, "P");
                }
            }
        } // end loop over wSet

        // BestWs
        int bestW_old = bfdOld.bestW;
        int bestW_new = bfdNew.bestW;
        
        // Build a 2×3 canvas
        TCanvas cCmp(Form("compareAll_s%d_ib%d", s,b), "", 2000,2400);
        cCmp.Divide(2,3,0.001, 0.001);

        // ---- Retrieve min/max of mgOld & mgNew => unify top plots
        auto [minOld, maxOld] = getMinMaxY(mgOld);
        auto [minNew, maxNew] = getMinMaxY(mgNew);
        double combinedMin = std::min(minOld, minNew);
        double combinedMax = std::max(maxOld, maxNew);
        // pad slightly
        combinedMin = std::max(0.0, 0.9*combinedMin);
        combinedMax *= 1.1;

        // 1) top-left => mgOld
        cCmp.cd(1);
        {
            int nGo = (mgOld->GetListOfGraphs() ? mgOld->GetListOfGraphs()->GetSize() : 0);
            if(nGo>0) {
                mgOld->Draw("A");
                mgOld->GetXaxis()->SetLimits(-2000,2000);
                mgOld->SetMinimum(combinedMin);
                mgOld->SetMaximum(combinedMax);
                gPad->SetGrid();
            } else {
                TLatex la; la.DrawLatexNDC(0.3,0.5,"No OLD data!");
            }
        }

        // 2) top-right => mgNew
        cCmp.cd(2);
        {
            int nGn = (mgNew->GetListOfGraphs() ? mgNew->GetListOfGraphs()->GetSize() : 0);
            if(nGn>0) {
                mgNew->Draw("A");
                mgNew->GetXaxis()->SetLimits(-2000,2000);
                mgNew->SetMinimum(combinedMin);
                mgNew->SetMaximum(combinedMax);
                gPad->SetGrid();
            } else {
                TLatex la; la.DrawLatexNDC(0.3,0.5,"No NEW data!");
            }
        }

        // 3) middle-left => ratio
        cCmp.cd(3);
        {
            gPad->SetGrid();
            int nr = (mgRatio->GetListOfGraphs() ?
                      mgRatio->GetListOfGraphs()->GetSize() : 0);
            if(nr>0) mgRatio->Draw("A");
            else {
                TLatex la; la.DrawLatexNDC(0.3,0.5,"No ratio data!");
            }
        }

        // 4) middle-right => build a legend with a simpler format:
        cCmp.cd(4);
        {
            // One approach: we create a "dummy" set of markers for each width
            // to show color, plus separate dummy markers for old/new shape.
            TLegend* leg = new TLegend(0.05,0.05,0.95,0.95);
            leg->SetNColumns(2);
            leg->SetBorderSize(1);
            leg->SetFillColorAlpha(kWhite,0.8);
            leg->SetTextSize(0.03);

            // for each w in wSet, add a single entry for color
            for (int w : wSet) {
                TMarker* mDummy = new TMarker(0,0,20); // shape doesn't matter
                mDummy->SetMarkerColor(colorForWidth(w));
                mDummy->SetMarkerStyle(20);
                leg->AddEntry(mDummy, Form("W=%d", w), "p");
            }

            // Then add separate markers to show "Old => circle(20), New => square(21), Ratio=>triangle(25)"
            TMarker *mOld = new TMarker(0,0,20);
            TMarker *mNew = new TMarker(0,0,21);
            TMarker *mRat = new TMarker(0,0,25);

            mOld->SetMarkerColor(kBlack);
            mNew->SetMarkerColor(kBlack);
            mRat->SetMarkerColor(kBlack);

            leg->AddEntry(mOld, "OLD data",   "p");
            leg->AddEntry(mNew, "NEW data",   "p");
            leg->AddEntry(mRat, "RATIO data", "p");

            leg->Draw();
        }

        // ----------------------------------------------------------------------
        // bottom-left => old bestW fit
        // ----------------------------------------------------------------------
        double kOld=0., oldA=0.;
        cCmp.cd(5);
        {
            auto itData = bfdOld.data.find(bestW_old);
            if (bestW_old<1 || itData==bfdOld.data.end() || itData->second.empty()) {
                TLatex la; la.DrawLatexNDC(0.3,0.5,"No OLD bestW data");
            }
            else
            {
                const auto &arr = itData->second;
                TGraphErrors* gE = new TGraphErrors((int)arr.size());

                double minY=1e9, maxY=-1e9;
                for(int i=0;i<(int)arr.size(); i++){
                    gE->SetPoint(i, arr[i].x, arr[i].y);
                    gE->SetPointError(i,0., arr[i].e);
                    if(arr[i].y<minY) minY=arr[i].y;
                    if(arr[i].y>maxY) maxY=arr[i].y;
                }
                double yLo= (minY>0 ? 0.8*minY : 0.);
                double yHi= 1.2*maxY;

                gPad->DrawFrame(-2000,yLo,2000,yHi,
                                Form("OLD bestW=%d;#Delta(mV);AvgADC", bestW_old));
                gE->SetMarkerStyle(20);
                gE->SetMarkerColor(kBlue+1);
                gE->Draw("P SAME");

                // re-draw exponential
                TF1* fExp = new TF1("fExp_old","[0]*exp([1]*x)",-1500,1500);
                fExp->SetLineColor(kRed);
                fExp->SetParameter(0, bfdOld.amplitude);
                fExp->SetParameter(1, bfdOld.kParam);
                fExp->Draw("SAME");

                oldA  = bfdOld.amplitude;
                kOld  = bfdOld.kParam;
                double eA = bfdOld.amplitudeErr;
                double eK = bfdOld.kParamErr;
                double resolutionK= (std::fabs(kOld)>1e-12 ?
                                     (eK/std::fabs(kOld))*100. : 0.);

                TLatex lat; lat.SetNDC(true);
                lat.SetTextSize(0.04);
                double yT=0.78;
                lat.DrawLatex(0.15,yT, Form("A=%.2f #pm %.2f", oldA,eA));
                yT -=0.06;
                lat.DrawLatex(0.15,yT, Form("k=%.4g #pm %.4g", kOld,eK));
                yT -=0.06;
                lat.DrawLatex(0.15,yT, Form("Res(k)=%.2f%%", resolutionK));
            }
        }

        // ----------------------------------------------------------------------
        // bottom-right => new bestW fit
        // ----------------------------------------------------------------------
        double kNew=0., newA=0.;
        cCmp.cd(6);
        {
            auto itData = bfdNew.data.find(bestW_new);
            if (bestW_new<1 || itData==bfdNew.data.end() || itData->second.empty()) {
                TLatex la; la.DrawLatexNDC(0.3,0.5,"No NEW bestW data");
            }
            else
            {
                const auto &arr = itData->second;
                TGraphErrors* gE = new TGraphErrors((int)arr.size());

                double minY=1e9, maxY=-1e9;
                for(int i=0;i<(int)arr.size(); i++){
                    gE->SetPoint(i, arr[i].x, arr[i].y);
                    gE->SetPointError(i,0., arr[i].e);
                    if(arr[i].y<minY) minY=arr[i].y;
                    if(arr[i].y>maxY) maxY=arr[i].y;
                }
                double yLo= (minY>0 ? 0.8*minY : 0.);
                double yHi= 1.2*maxY;

                gPad->DrawFrame(-2000,yLo,2000,yHi,
                                Form("NEW bestW=%d;#Delta(mV);AvgADC", bestW_new));
                gE->SetMarkerStyle(21);
                gE->SetMarkerColor(kRed+1);
                gE->Draw("P SAME");

                TF1* fExp = new TF1("fExp_new","[0]*exp([1]*x)",-1500,1500);
                fExp->SetLineColor(kBlue);
                fExp->SetParameter(0, bfdNew.amplitude);
                fExp->SetParameter(1, bfdNew.kParam);
                fExp->Draw("SAME");

                newA = bfdNew.amplitude;
                kNew = bfdNew.kParam;
                double eA = bfdNew.amplitudeErr;
                double eK = bfdNew.kParamErr;
                double resolutionK= (std::fabs(kNew)>1e-12 ?
                                     (eK/std::fabs(kNew))*100. : 0.);

                TLatex lat; lat.SetNDC(true);
                lat.SetTextSize(0.04);
                double yT=0.78;
                lat.DrawLatex(0.15,yT, Form("A=%.2f #pm %.2f", newA,eA));
                yT -=0.06;
                lat.DrawLatex(0.15,yT, Form("k=%.4g #pm %.4g", kNew,eK));
                yT -=0.06;
                lat.DrawLatex(0.15,yT, Form("Res(k)=%.2f%%", resolutionK));
            }
        }

        // Store final (k) + amplitude in your existing arrays
        kRatios.push_back({ s,b, kOld, kNew });
        amplitudeRatios.push_back({ s,b, oldA, newA });

        // Save final multiW overlay
        const std::string saveDir =
            "/Users/patsfan753/Desktop/PMTgainsAna/output/PerIBOverlays/compare";
        gSystem->Exec(Form("mkdir -p %s", saveDir.c_str()));
        
        std::string outName = Form("%s/compareAll_s%d_ib%d.png",
                                   saveDir.c_str(), s, b);
        cCmp.SaveAs(outName.c_str());
    } // end loop over oldFitMap
    
    std::cout << "[INFO] Created " << validIB
              << " IB multiW overlay plots in " << multiW << "\n\n";
    
    // Additional summary
    {
        double sumKold=0., sumKnew=0.;
        double minKold=1e9, maxKold=-1e9;
        double minKnew=1e9, maxKnew=-1e9;
        int count=0;
        
        for (auto &kk : kRatios) {
            if (kk.kOld>0) {
                sumKold += kk.kOld;
                if (kk.kOld<minKold) minKold= kk.kOld;
                if (kk.kOld>maxKold) maxKold= kk.kOld;
            }
            if (kk.kNew>0) {
                sumKnew += kk.kNew;
                if (kk.kNew<minKnew) minKnew= kk.kNew;
                if (kk.kNew>maxKnew) maxKnew= kk.kNew;
            }
            count++;
        }
        double avgOld = (count>0 ? sumKold/count : 0.);
        double avgNew = (count>0 ? sumKnew/count : 0.);
        
        std::cout << "[DEBUG] Additional summary of k values:\n"
                  << "  # IB = " << count << "\n"
                  << "  OLD: avg(k) = " << avgOld
                  << ", min = " << minKold
                  << ", max = " << maxKold << "\n"
                  << "  NEW: avg(k) = " << avgNew
                  << ", min = " << minKnew
                  << ", max = " << maxKnew << "\n\n";
    }
}


/*****************************************************************************
 *  Low-level helper : draws every (sector,ib) overlay for ONE data-set
 *  - The color of each curve depends on the test-pulse width `w`
 *  - We only draw the exponential fit (in red) for the bestW, using the
 *    amplitude/kParam from the stored BoardFitData.
 *  - We also produce a multi-page "tabulated best-fits" (6×5) summary, with
 *    each board's bestW curve + exponential fit in a sub-pad.
 *****************************************************************************/

static void produceOneSetOfOverlays(
    const std::map< std::pair<int,int>, BoardFitData > &fitMap,
    const std::string &tag,                       // "OLD" / "NEW"
    const std::string &outDirBase,                // e.g. …/overlays/old or …/overlays/new
    const std::map<int,Color_t> &widthColour      // fixed colour table for widths
)
{
    /*--------------------------------------------------------------------*/
    /* 1) ensure the output directory exists                              */
    /*--------------------------------------------------------------------*/
    gSystem->Exec(Form("mkdir -p %s", outDirBase.c_str()));

    /*--------------------------------------------------------------------*/
    /* 2) track boards with amplitude>6000, for final summary             */
    /*--------------------------------------------------------------------*/
    struct HighAmpInfo {
        int s;
        int ib;
        double A;
        double eA;
    };
    std::vector<HighAmpInfo> highAmpList;

    /*--------------------------------------------------------------------*/
    /* We'll also collect a struct with each board's "bestW" data so we   */
    /* can produce the 6×5 tabulated best-fits at the end.                */
    /*--------------------------------------------------------------------*/
    struct BestFitBoard {
        int s, ib;
        int bestW;
        double amplitude, amplitudeErr;
        double kParam, kParamErr;
        // We'll also store the data vector for bestW
        std::vector<XYerr> bestWdata;
    };
    std::vector<BestFitBoard> tabulated; // we fill this as we go

    /*--------------------------------------------------------------------*/
    /* 3) loop over every (sector,ib) in fitMap => produce normal overlays*/
    /*--------------------------------------------------------------------*/
    for (const auto &pr : fitMap)
    {
        int s  = pr.first.first;
        int ib = pr.first.second;

        // optional skip if "bad board"
        if (isBadBoard(s,ib))
            continue;

        const BoardFitData &bfd = pr.second;
        if (bfd.data.empty())
            continue;  // no data => skip

        // The best test-pulse width from BoardFitData
        int bestW = bfd.bestW;
        if (bestW <= 0) {
            // if bestW isn't valid, skip
            continue;
        }

        /*----------------------------------------------------------------*/
        /* 4) Build a TCanvas + TMultiGraph for all widths on that IB     */
        /*----------------------------------------------------------------*/
        TCanvas c(Form("%sGain_S%d_IB%d", tag.c_str(), s, ib),
                  Form("%sGain_S%d_IB%d", tag.c_str(), s, ib),
                  1200,900);
        c.SetRightMargin(0.28);
        c.SetLeftMargin(0.15);

        TMultiGraph *mg = new TMultiGraph();
        mg->SetTitle(Form("Gain curves (%s), Sector %d, IB %d;#Delta(mV);Avg ADC",
                          tag.c_str(), s, ib));

        // A legend on the right side
        TLegend leg(0.72,0.10,0.98,0.90);
        leg.SetNColumns(2);
        leg.SetTextSize(0.03);
        leg.SetBorderSize(1);
        leg.SetFillColorAlpha(kWhite,0.7);

        TGraphErrors *grBest = nullptr; // We'll store the TGraph for bestW

        /*----------------------------------------------------------------*/
        /* 5) For each width => build TGraphErrors => add to mg           */
        /*----------------------------------------------------------------*/
        for (auto &kv : bfd.data)
        {
            int w = kv.first;
            const auto &xyVec = kv.second;
            if (xyVec.empty())
                continue;

            // create TGraphErrors
            TGraphErrors *gr = new TGraphErrors((int)xyVec.size());
            for (int i=0; i<(int)xyVec.size(); i++){
                gr->SetPoint(i, xyVec[i].x, xyVec[i].y);
                gr->SetPointError(i, 0., xyVec[i].e);
            }

            // pick color from widthColour map if found; else fallback
            auto itCol = widthColour.find(w);
            Color_t col = (itCol != widthColour.end() ? itCol->second : kBlack);

            // highlight bestW in black
            if (w == bestW) {
                gr->SetMarkerStyle(21);
                gr->SetMarkerColor(kBlack);
                gr->SetLineColor(kBlack);
                grBest = gr;
            }
            else {
                // normal widths => marker style=20, color=col
                gr->SetMarkerStyle(20);
                gr->SetMarkerColor(col);
                gr->SetLineColor(col);
            }

            gr->SetNameTitle(Form("w%d", w), Form("W=%d ns", w));
            mg->Add(gr, "P");
            leg.AddEntry(gr, Form("W=%d", w), "lp");
        }

        /*----------------------------------------------------------------*/
        /* 6) draw the TMultiGraph, set y-range, etc.                     */
        /*----------------------------------------------------------------*/
        mg->Draw("A");
        mg->GetXaxis()->SetLimits(-2000,2000);
        mg->SetMinimum(0.);
        gPad->SetGrid();

        // expand yMax by ~20%
        double yMax = (mg->GetHistogram() ? mg->GetHistogram()->GetMaximum() : 0.);
        mg->SetMaximum(1.2 * yMax);

        leg.Draw();

        /*----------------------------------------------------------------*/
        /* 7) If we have bestW => draw the exponential fit in red         */
        /*----------------------------------------------------------------*/
        if (grBest) {
            // We'll use the amplitude & slope from bfd
            // We'll define a TF1 with that form & draw it
            TF1 fExp("fExp","[0]*exp([1]*x)",-1500,1500);
            fExp.SetLineColor(kRed);

            fExp.SetParameter(0, bfd.amplitude);
            fExp.SetParameter(1, bfd.kParam);

            double parErr[2] = { bfd.amplitudeErr, bfd.kParamErr };
            fExp.SetParErrors(parErr);

            fExp.Draw("SAME");

            // Add textual annotations near the top-left
            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextSize(0.028);
            double y0 = 0.85;

            lat.DrawLatex(0.16, y0, Form("BestW = %d ns", bestW));
            y0 -= 0.05;
            lat.DrawLatex(0.16, y0,
                          Form("A = %.1f #pm %.1f", bfd.amplitude, bfd.amplitudeErr));
            y0 -= 0.05;
            lat.DrawLatex(0.16, y0,
                          Form("k = %.3g #pm %.3g", bfd.kParam, bfd.kParamErr));

            // track high amplitude boards
            if (bfd.amplitude > 6000.) {
                HighAmpInfo hi { s, ib, bfd.amplitude, bfd.amplitudeErr };
                highAmpList.push_back(hi);
            }
        }

        /*----------------------------------------------------------------*/
        /* 8) Save as .png                                                */
        /*----------------------------------------------------------------*/
        std::string outPath = Form("%s/S%d_IB%d.png", outDirBase.c_str(), s, ib);
        c.SaveAs(outPath.c_str());

        // clean up
        delete mg;

        /*----------------------------------------------------------------*/
        /* 9) Also store info for the “tabulated best-fits” 6×5           */
        /*----------------------------------------------------------------*/
        auto itBW = bfd.data.find(bestW);
        if (itBW != bfd.data.end()) {
            BestFitBoard bfb;
            bfb.s = s;
            bfb.ib = ib;
            bfb.bestW = bestW;
            bfb.amplitude     = bfd.amplitude;
            bfb.amplitudeErr  = bfd.amplitudeErr;
            bfb.kParam        = bfd.kParam;
            bfb.kParamErr     = bfd.kParamErr;
            bfb.bestWdata     = itBW->second;  // the XYerr vector

            tabulated.push_back(bfb);
        }
    } // end board loop

    /*--------------------------------------------------------------------*/
    /* 10) final summary of boards with amplitude>6000                    */
    /*--------------------------------------------------------------------*/
    std::cout << "\n\033[1;33m==========  SUMMARY OF " << tag
              << " FIT AMPLITUDES >6000  ==========\033[0m\n";
    if (highAmpList.empty()) {
        std::cout << "\033[1;32mNo boards had amplitude above 6000.\033[0m\n";
    }
    else {
        std::cout << "\033[1;36mFound " << highAmpList.size()
                  << " board(s) with amplitude>6000.\033[0m\n"
                  << "   Sector | IB  |  Fitted A ± eA\n"
                  << "   -------+-----+------------------\n";
        for (auto &hi : highAmpList) {
            std::cout << "   " << std::setw(6) << hi.s << " | "
                      << std::setw(3) << hi.ib << " | "
                      << Form("%9.1f ± %4.1f", hi.A, hi.eA)
                      << "\n";
        }
    }
    std::cout << "\033[1;33m==========================================================\033[0m\n"
              << "[INFO] → per-IB overlays for " << tag
              << " saved in " << outDirBase << '\n';

    /*--------------------------------------------------------------------*/
    /* 11) 6×5 tabulated best-fits for each IB                            */
    /*--------------------------------------------------------------------*/
    std::string tableFolder;
    if (tag=="OLD")
        tableFolder = "tabulatedBestFitsOld";
    else
        tableFolder = "tabulatedBestFitsNew";

    std::string outFolder = Form("%s/%s", outDirBase.c_str(), tableFolder.c_str());
    gSystem->Exec(Form("mkdir -p %s", outFolder.c_str()));

    const int boardsPerCanvas = 6*5;  // 30 subpads
    int totalBoards = (int)tabulated.size();
    int pageIndex=0;

    for(int i=0; i < totalBoards; i += boardsPerCanvas)
    {
        int chunkSize = std::min(boardsPerCanvas, totalBoards - i);

        TCanvas cTab(Form("tableOfBestFits_%s_page%d", tag.c_str(), pageIndex),
                     Form("BestFits %s – page #%d", tag.c_str(), pageIndex),
                     2000, 1800);
        cTab.Divide(6,5, 0.001, 0.001);  // 6 across, 5 down

        for(int j=0; j<chunkSize; j++){
            int padIndex = j+1;
            cTab.cd(padIndex);

            const BestFitBoard &bfb = tabulated[i + j];

            // draw a mini-plot for the "bestW" only
            TGraphErrors *gE = new TGraphErrors((int)bfb.bestWdata.size());
            double minY=1e9, maxY=-1e9;
            for(int k=0; k<(int)bfb.bestWdata.size(); k++){
                double xx = bfb.bestWdata[k].x;
                double yy = bfb.bestWdata[k].y;
                double ee = bfb.bestWdata[k].e;
                gE->SetPoint(k,xx,yy);
                gE->SetPointError(k,0.,ee);
                if(yy<minY) minY=yy;
                if(yy>maxY) maxY=yy;
            }

            double yLo= (minY>0 ? 0.8*minY : 0.);
            double yHi= 1.2*maxY;

            gPad->DrawFrame(-2000,yLo,2000,yHi,
                            Form("S%d.IB%d => bestW=%d;#Delta(mV);ADC",
                                 bfb.s, bfb.ib, bfb.bestW));
            gE->SetMarkerStyle(21);
            gE->SetMarkerColor(kBlack);
            gE->Draw("P SAME");

            // draw the best-fit line from the stored amplitude/k
            TF1 fExp("fExpBW","[0]*exp([1]*x)",-1500,1500);
            fExp.SetLineColor(kRed);
            fExp.SetParameter(0, bfb.amplitude);
            fExp.SetParameter(1, bfb.kParam);

            double parErr[2] = { bfb.amplitudeErr, bfb.kParamErr };
            fExp.SetParErrors(parErr);

            fExp.Draw("SAME");

            // small annotation
            {
                TLatex ltx;
                ltx.SetNDC(true);
                ltx.SetTextSize(0.04);
                double yAnn=0.75;
                ltx.DrawLatex(0.15,yAnn, Form("A=%.1f #pm %.1f", bfb.amplitude, bfb.amplitudeErr));
                yAnn -=0.06;
                ltx.DrawLatex(0.15,yAnn, Form("k=%.3g #pm %.3g", bfb.kParam, bfb.kParamErr));
            }
        }
        std::string outImg = Form("%s/bestFitsPage_%d.png", outFolder.c_str(), pageIndex);
        cTab.SaveAs(outImg.c_str());
        pageIndex++;
    }

    std::cout << "[INFO] => Created " << pageIndex
              << " tabulated best-fit pages for " << tag
              << " in " << outFolder << std::endl;
}


/*****************************************************************************
 *   High-level wrapper – decides a fixed colour per width *once*, then
 *   calls the low-level routine for OLD and NEW so colours match.
 *****************************************************************************/
void producePerIBOverlays(
    const std::map< std::pair<int,int>, BoardFitData > &oldFitMap,
    const std::map< std::pair<int,int>, BoardFitData > &newFitMap,
    const std::string &baseOutDir)
{
    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));

    /*--------------------------------------------------------------------*/
    /* 1) gather all widths from old & new to unify color assignment      */
    /*--------------------------------------------------------------------*/
    std::set<int> allW;
    auto collectW = [&allW](const auto &m){
        for (const auto &pr : m) {
            for (const auto &kv : pr.second.data) {
                allW.insert(kv.first);
            }
        }
    };
    collectW(oldFitMap);
    collectW(newFitMap);

    /*--------------------------------------------------------------------*/
    /* 2) build deterministic (width -> colour) mapping                   */
    /*--------------------------------------------------------------------*/
    std::map<int,Color_t> widthColour;
    size_t idx=0;
    for (int w : allW) {
        widthColour[w] = kPalette[idx % kPalette.size()];
        ++idx;
    }

    /*--------------------------------------------------------------------*/
    /* 3) Make subfolders and delegate to the drawing routine             */
    /*--------------------------------------------------------------------*/
    const std::string outOld = baseOutDir + "/old";
    const std::string outNew = baseOutDir + "/new";

    produceOneSetOfOverlays(oldFitMap, "OLD", outOld, widthColour);
    produceOneSetOfOverlays(newFitMap, "NEW", outNew, widthColour);
}



void doAnalysis(
    const std::string &oldRunCSV,
    const std::string &oldInputDir,
    const std::string &newRunCSV,
    const std::string &newInputDir,
    const std::string &vopCSV
)
{
    //----------------------------------------------------------------------
    // 0) Basic prints + ensure output dirs
    //----------------------------------------------------------------------
    std::cout << "\n[DEBUG] CompareOldVsNewData() => self-contained\n"
              << " oldRunCSV="   << oldRunCSV    << "\n"
              << " oldInputDir=" << oldInputDir  << "\n"
              << " newRunCSV="   << newRunCSV    << "\n"
              << " newInputDir=" << newInputDir  << "\n"
              << " vopCSV="      << vopCSV       << "\n\n";


    //----------------------------------------------------------------------
    // 1) Build best-pulse-width fits (both old & new) => oldFitMap, newFitMap
    //    using the function you specified:
    //----------------------------------------------------------------------
    std::map< std::pair<int,int>, BoardFitData > oldFitMap;
    std::map< std::pair<int,int>, BoardFitData > newFitMap;

    buildBestPulseWidthMap(
        oldRunCSV,
        oldInputDir,
        newRunCSV,
        newInputDir,
        vopCSV,
        oldFitMap,
        newFitMap
    );
    
    //----------------------------------------------------------------------
    // 2) Per-board multi-width overlays  (NEW!)
    //----------------------------------------------------------------------
    const std::string overlayDir =
        "/Users/patsfan753/Desktop/PMTgainsAna/output/PerIBOverlays/compare";

    // create the directory once
    gSystem->Exec(Form("mkdir -p %s", overlayDir.c_str()));

    std::vector<Kval>   kRatios;         // filled inside the routine
    std::vector<AmpVal> amplitudeRatios; //     »       »       »

    produceMultiWOverlayPlots(oldFitMap,
                              newFitMap,
                              /*multiW =*/overlayDir,
                              /*summary=*/overlayDir,   // only used for a log line
                              kRatios,
                              amplitudeRatios);

    produceCompare32vsBest(oldFitMap, newFitMap);

    //--------------------------------------------------------------
    produceColorCodedMapsAndExtra(
        oldFitMap,
        newFitMap,
        "/Users/patsfan753/Dekstop/PMTgainsAna/output/2Dmaps"
    );


    auto kStats = produceKAndAmplitudePlots(oldFitMap, newFitMap);
    double muK_new   = kStats.first;
    double sigmaK_new= kStats.second;
    

    const std::string summaryDir =
        "/Users/patsfan753/Desktop/PMTgainsAna/output/summary_ampAndk_vals";
    const std::string ampPNG =
        summaryDir + "/amplitudeCheck_outliers_twopanels.png";

    outOfSigmaAmplitudeCheck(newFitMap,
                             oldFitMap,
                             muK_new,
                             sigmaK_new,
                             ampPNG);


    
    producePerIBOverlays(
        oldFitMap,
        newFitMap,
        "/Users/patsfan753/Desktop/PMTgainsAna/output/PerIBOverlays"  // e.g. "/Users/.../PerIBOverlays"
    );

}




//--------------------------------------------------------------------------------
// main analysis
//--------------------------------------------------------------------------------
void PMTanalysisTestPulse()
{
    std::string basePath   = "/Users/patsfan753/Desktop/PMTgainsAna/data/";
    
    doAnalysis(
        basePath + "04-11-25-run-info.csv",
        basePath + "input",
        basePath + "04-26-25-run-info.csv",
        basePath + "inputNew",
        basePath + "vop.csv"
    );
 
  // End with the original final message:
  std::cout<<"[INFO] ---------- ANALYSIS COMPLETE ---------.\n";
}
