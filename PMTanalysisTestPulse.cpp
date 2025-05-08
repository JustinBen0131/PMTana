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

// Forward-declare the XYerr struct if needed:
struct XYerr {
    double x;
    double y;
    double e;
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


// ----------------------------------------------------------------------
// buildBestPulseWidthMap(...) => fills oldFitMap, newFitMap
// ----------------------------------------------------------------------
void buildBestPulseWidthMap(const std::string &oldRunCSV,
                            const std::string &oldInputDir,
                            const std::string &newRunCSV,
                            const std::string &newInputDir,
                            const std::string &vopCSV,
                            // OUT references => your final maps
                            std::map< std::pair<int,int>, BoardFitData > &oldFitMap,
                            std::map< std::pair<int,int>, BoardFitData > &newFitMap)
{
    // 1) read old/new CSV => build runDB
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
    auto isOld = [&](int runID){ return (oldSet.count(runID)>0); };

    // 2) read VOP => build (sector,ib) map
    std::vector<TowerDefault> allTows = readVopCSV(vopCSV);
    std::map<std::pair<int,int>, bool> sectorIBmap;
    for(auto &td : allTows){
        sectorIBmap[{td.sector, td.ib}] = true;
    }

    // 3) We'll store data in dataMap[0 or 1][(sector,ib)][width]
    std::map<int, std::map< std::pair<int,int>, std::map<int, std::vector<XYerr>> >> dataMap;

    // Helper => fill data from a single run
    auto processRun = [&](int runID, const std::string &dir, bool oldNotNew){
        auto it = runDB.find(runID);
        if(it==runDB.end()) return;

        RunInfo &ri = it->second;
        TH2F* h2 = getH2CEMC(dir, runID);
        if(!h2) return;
        zeroOutBadBoards(h2);

        for(auto &sb: sectorIBmap) {
            int s= sb.first.first;
            int b= sb.first.second;
            if(isBadBoard(s,b)) continue;

            // find default offset
            double defOff=0.;
            {
                double sum=0.;
                int ct=0;
                for(auto &td : allTows){
                    if(td.sector==s && td.ib==b){
                        sum+= td.defaultOffset;
                        ct++;
                    }
                }
                defOff = (ct>0? sum/ct : 0.);
            }
            double delta=0.;
            if(!ri.isOriginal){
                delta = ri.numericOffset - defOff;
            }
            // else delta=0 => means "original"

            // gather average ADC
            double sumVal=0., sumErr2=0.;
            int count=0;
            for(auto &td: allTows){
                if(td.sector==s && td.ib==b) {
                    int bx= h2->GetXaxis()->FindBin(td.iphi);
                    int by= h2->GetYaxis()->FindBin(td.ieta);
                    double val = h2->GetBinContent(bx,by);
                    double err = h2->GetBinError(bx,by);
                    if(val<0) continue;
                    sumVal+= val;
                    sumErr2+=(err*err);
                    count++;
                }
            }
            if(count<1) continue;

            double meanADC = sumVal / count;
            double meanErr = std::sqrt(sumErr2)/ count;
            if(meanADC<=0) continue;

            XYerr ex{ delta, meanADC, meanErr };
            int idx= (oldNotNew? 0 : 1);
            dataMap[idx][{s,b}][ ri.width ].push_back(ex);
        }
        delete h2;
    };

    // Process old + new runs
    for(auto &rr : oldRuns) {
        processRun(rr.runNumber, oldInputDir, true);
    }
    for(auto &rr : newRuns){
        processRun(rr.runNumber, newInputDir, false);
    }

    // 4) find bestW for OLD => offset=0 near 3000
    std::map<std::pair<int,int>,int> bestWmap;
    {
        // offset=0 => store ADC
        std::map< std::pair<int,int>, std::map<int,std::vector<double>> > zeroStore;
        for(auto &sbItem : dataMap[0]){
            auto &Wmap = sbItem.second;
            for(auto &wItem : Wmap){
                int w = wItem.first;
                for(auto &xy: wItem.second){
                    if(std::fabs(xy.x)<1e-9 && xy.y>0){
                        zeroStore[sbItem.first][w].push_back(xy.y);
                    }
                }
            }
        }
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
                double mean= sum/vals.size();
                double df= std::fabs(mean-3000.);
                if(df< bestDiff){
                    bestDiff= df;
                    bestW= w;
                }
            }
            bestWmap[{s,b}] = bestW;
        }
    }

    // 5) bestW for NEW => offset=0 near 3000
    std::map<std::pair<int,int>,int> bestWmapNew;
    {
        std::map< std::pair<int,int>, std::map<int,std::vector<double>> > zeroNew;
        for(auto &sbN: dataMap[1]){
            auto &Wmap= sbN.second;
            for(auto &kv: Wmap){
                int w= kv.first;
                for(auto &xy : kv.second){
                    if(std::fabs(xy.x)<1e-9 && xy.y>0){
                        zeroNew[sbN.first][w].push_back(xy.y);
                    }
                }
            }
        }
        for(auto &xx : zeroNew){
            int s= xx.first.first;
            int b= xx.first.second;
            double bestDiff=1e9;
            int bestW=-1;
            for(auto &wvals : xx.second){
                if(wvals.second.empty()) continue;
                double sum=0.;
                for(double val: wvals.second) sum+= val;
                double mean= sum / wvals.second.size();
                double df= std::fabs(mean-3000.);
                if(df< bestDiff){
                    bestDiff= df;
                    bestW   = wvals.first;
                }
            }
            bestWmapNew[{s,b}] = bestW;
        }
    }

    // 6) Exponential fit => store in oldFitMap / newFitMap
    auto fitBestWidth = [&](bool oldNotNew, int sector, int ib, int bestW){
        BoardFitData result;
        result.isOld= oldNotNew;
        result.sector= sector;
        result.ib= ib;
        result.bestW= bestW;

        int idx= (oldNotNew? 0:1);
        auto &arr = dataMap[idx][{sector, ib}][bestW];
        if(arr.empty()) return result;

        // copy all widths => so we have them in .data
        for(auto &kv: dataMap[idx][{sector, ib}]){
            result.data[kv.first] = kv.second;
        }

        // build TGraphErrors
        TGraphErrors *gE = new TGraphErrors((int)arr.size());
        for(int i=0; i<(int)arr.size(); i++){
            gE->SetPoint(i, arr[i].x, arr[i].y);
            gE->SetPointError(i, 0., arr[i].e);
        }
        TF1 fExp("fExp","[0]*exp([1]*x)", -1500,1500);
        double x0,y0;
        gE->GetPoint(0,x0,y0);
        double yFirst= (y0>0.? y0:1000.);
        fExp.SetParameters(yFirst,0.0005);

        TFitResultPtr fr = gE->Fit(&fExp,"S E R Q");

        double A  = fExp.GetParameter(0);
        double eA = fExp.GetParError(0);
        double k  = fExp.GetParameter(1);
        double eK = fExp.GetParError(1);

        result.amplitude    = A;
        result.amplitudeErr = eA;
        result.kParam       = k;
        result.kParamErr    = eK;

        delete gE;
        return result;
    };

    // For old => fill oldFitMap
    for(auto &item: dataMap[0]){
        int s= item.first.first;
        int b= item.first.second;
        int w= bestWmap[{s,b}];
        BoardFitData bfd = fitBestWidth(true, s,b,w);
        oldFitMap[{s,b}] = bfd;
    }
    // For new => fill newFitMap
    for(auto &item: dataMap[1]){
        int s= item.first.first;
        int b= item.first.second;
        int w= bestWmapNew[{s,b}];
        BoardFitData bfd= fitBestWidth(false, s,b,w);
        newFitMap[{s,b}] = bfd;
    }
}


void produceNewPerIBOverlays(
    const std::map< std::pair<int,int>, BoardFitData > & newFitMap,
    const std::vector< /* your amplitude struct */ AmpVal > & amplitudeRatios,
    const std::map<std::pair<int,int>, bool> & sectorIBmap,
    const std::string &outDir
)
{
    // 1) Ensure the output folder
    gSystem->Exec(Form("mkdir -p %s", outDir.c_str()));

    // A small struct to track boards with amplitude>6000
    struct HighAmpInfo {
        int sector;
        int ib;
        double A_fitted;    // from newFitMap
        double eA_fitted;   // from newFitMap
        double A_stored;    // from amplitudeRatios
    };
    std::vector<HighAmpInfo> highAmpList;

    // We'll define a simple color palette for each width. You can reuse your existing:
    std::vector<Color_t> colorVec = {
        kBlue, kGreen+2, kMagenta, kCyan+2,
        kOrange+1, kViolet, kGray+2, kBlack
    };

    // 2) Loop over each (sector,IB) in sectorIBmap
    for (auto &sb : sectorIBmap)
    {
        int s = sb.first.first;
        int b = sb.first.second;
        if (isBadBoard(s,b)) continue; // skip known bad boards

        // Check if this board is in newFitMap
        auto itMap = newFitMap.find({s,b});
        if (itMap == newFitMap.end()) continue; // no data for that board

        const BoardFitData &bfd = itMap->second;
        // bfd.data[w] => your XY points
        // bfd.bestW   => best width
        // bfd.amplitude, bfd.amplitudeErr, bfd.kParam, bfd.kParamErr => final fit

        // If there's nothing in bfd.data, skip
        if (bfd.data.empty()) continue;

        // 3) Create a big TCanvas for clarity
        TCanvas cGain(Form("NewGain_S%d_IB%d", s,b),
                      Form("NewGain_S%d_IB%d", s,b),
                      1200, 900);
        cGain.SetRightMargin(0.28);
        cGain.SetLeftMargin(0.15);

        // A TMultiGraph to overlay all widths
        TMultiGraph *mg = new TMultiGraph();
        mg->SetTitle(Form("Gain Curves (NEW) for Sector %d, IB %d;#Delta(mV);Avg ADC", s,b));

        // Build a TLegend on the right side
        TLegend *leg = new TLegend(0.72, 0.10, 0.98, 0.90);
        leg->SetNColumns(2);
        leg->SetTextSize(0.03);
        leg->SetBorderSize(1);
        leg->SetFillColorAlpha(kWhite, 0.7);

        // We'll track the TGraph of the bestW (to draw an exponential with it)
        TGraphErrors *grBest = nullptr;

        // 4) For each width in bfd.data => build a TGraphErrors
        int colorIdx=0;
        for (auto &kv : bfd.data) {
            int w = kv.first;
            const auto &xyVec = kv.second;
            if (xyVec.empty()) continue;

            // build TGraphErrors
            TGraphErrors *gr = new TGraphErrors((int)xyVec.size());
            for(int i=0; i<(int)xyVec.size(); i++){
                gr->SetPoint(i, xyVec[i].x, xyVec[i].y);
                gr->SetPointError(i, 0., xyVec[i].e);
            }

            // decide color
            Color_t col = colorVec[colorIdx % colorVec.size()];
            colorIdx++;

            // highlight the bestW in black
            if (w == bfd.bestW) {
                gr->SetMarkerStyle(21);
                gr->SetMarkerColor(kBlack);
                gr->SetLineColor(kBlack);
                grBest = gr;
            } else {
                gr->SetMarkerStyle(20);
                gr->SetMarkerColor(col);
                gr->SetLineColor(col);
            }

            // add to legend
            gr->SetNameTitle(Form("gW%d", w), Form("W=%d ns", w));
            leg->AddEntry(gr, Form("W=%d", w), "lp");
            mg->Add(gr, "P");
        }

        // 5) Draw the multigraph
        mg->Draw("A");
        mg->GetXaxis()->SetLimits(-2000, 2000);
        mg->SetMinimum(0.);
        cGain.Update();

        double curMax = 0.;
        if(mg->GetHistogram()) {
            curMax = mg->GetHistogram()->GetMaximum();
        }
        mg->SetMaximum(1.2 * curMax);

        gPad->SetGrid();
        leg->Draw();

        // 6) Instead of re-fitting, we define a TF1 with the *stored* amplitude/k
        if (grBest) {
            TF1 *fExp = new TF1("fExp","[0]*exp([1]*x)", -1500,1500);
            fExp->SetLineColor(kRed);
            fExp->SetParameter(0, bfd.amplitude); // from newFitMap
            fExp->SetParError(0, bfd.amplitudeErr);
            fExp->SetParameter(1, bfd.kParam);
            fExp->SetParError(1, bfd.kParamErr);

            // optionally just "draw" it
            fExp->Draw("SAME");

            // We'll also place a few lines of text
            {
                TLatex lat;
                lat.SetNDC(true);
                lat.SetTextSize(0.028);
                double yLine = 0.85;

                lat.DrawLatex(0.16,yLine, Form("BestW = %d ns", bfd.bestW));
                yLine -= 0.04;
                lat.DrawLatex(0.16,yLine, Form("A=%.1f #pm %.1f", bfd.amplitude, bfd.amplitudeErr));
                yLine -= 0.04;
                lat.DrawLatex(0.16,yLine, Form("k=%.3g #pm %.3g", bfd.kParam, bfd.kParamErr));
            }

            // Also see if amplitude>6000 => track in highAmpList
            if (bfd.amplitude > 6000.) {
                // find the stored amplitude from amplitudeRatios
                double storedA_new=0.;
                for (auto &av : amplitudeRatios) {
                    if(av.sector==s && av.ib==b) {
                        storedA_new = av.ampNew;
                        break;
                    }
                }
                HighAmpInfo hi{ s, b, bfd.amplitude, bfd.amplitudeErr, storedA_new };
                highAmpList.push_back(hi);
            }
        }

        // 7) Save the canvas
        std::string outName = Form("%s/S%d_IB%d.png", outDir.c_str(), s,b);
        cGain.SaveAs(outName.c_str());
        delete mg; // cleanup
    }

    // 8) Print a summary of boards that had amplitude>6000
    std::cout << "\n\033[1;33m==========  SUMMARY OF NEW FIT AMPLITUDES >6000  ==========\033[0m\n";
    if (highAmpList.empty()) {
        std::cout << "\033[1;32mNo boards had amplitude above 6000.\033[0m\n";
    } else {
        std::cout << "\033[1;36mFound " << highAmpList.size()
                  << " board(s) with amplitude>6000.\033[0m\n"
                  << "   Sector | IB  |  Fitted A ± eA   |  A from amplitudeRatios\n"
                  << "   -------+-----+------------------+------------------------\n";
        for (auto &hi : highAmpList) {
            std::cout << "   "
                      << std::setw(6) << hi.sector << " | "
                      << std::setw(3) << hi.ib     << " | "
                      << Form("%9.1f ± %4.1f", hi.A_fitted, hi.eA_fitted)
                      << " | "
                      << Form("%9.1f  (stored)", hi.A_stored)
                      << "\n";
        }
    }
    std::cout << "\033[1;33m==========================================================\033[0m\n\n";
    std::cout << "[INFO] => Created per-IB overlays in " << outDir << "\n";
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
    // Just like your original code...
    std::cout << "[DEBUG] Step5 => multi-W overlay plots...\n";
    
    // We track how many IBs actually produce a multiW plot:
    int validIB=0;
    
    // We pick which widths to overlay. This matches your original {26..32} set.
    std::set<int> wSet{26,27,28,29,30,31,32};
    
    // Hard-coded colors for each overlay
    Color_t cols[8] = {kBlue, kRed, kGreen+2, kMagenta,
                       kCyan+2, kOrange+7, kBlack, kGray+2};
    
    // Loop over all boards in oldFitMap (or newFitMap).
    // If a board is in oldFitMap but not in newFitMap, presumably skip.
    for (auto &oldPair : oldFitMap)
    {
        int s = oldPair.first.first;
        int b = oldPair.first.second;
        
        // Check if newFitMap also has that board
        auto itNew = newFitMap.find({s,b});
        if (itNew == newFitMap.end()) {
            // Board not found in newFitMap => skip
            continue;
        }
        
        // Extract the BoardFitData
        const BoardFitData &bfdOld = oldPair.second;  // old
        const BoardFitData &bfdNew = itNew->second;   // new
        
        // If either .data is empty, skip
        if (bfdOld.data.empty() && bfdNew.data.empty()) continue;
        
        validIB++;
        
        // Build TMultiGraphs: mgOld, mgNew, mgRatio
        TMultiGraph* mgOld   = new TMultiGraph();
        mgOld->SetTitle(Form("OLD (S%d,IB%d);#Delta(mV);AvgADC", s,b));
        
        TMultiGraph* mgNew   = new TMultiGraph();
        mgNew->SetTitle(Form("NEW (S%d,IB%d);#Delta(mV);AvgADC", s,b));
        
        TMultiGraph* mgRatio = new TMultiGraph();
        mgRatio->SetTitle("New/Old;#Delta(mV);Ratio");
        
        // We'll gather each width from wSet => plot old vs. new
        int idx=0;
        for (int w : wSet)
        {
            // For OLD => bfdOld.data[w], for NEW => bfdNew.data[w]
            auto itOldW = bfdOld.data.find(w);
            auto itNewW = bfdNew.data.find(w);
            
            // If old does not have that width, skip
            // If new does not have that width, skip
            // (Or you can do partial skipping logic if you prefer.)
            const std::vector<XYerr>* vecO = nullptr;
            const std::vector<XYerr>* vecN = nullptr;
            
            if (itOldW != bfdOld.data.end())
                vecO = &itOldW->second; // pointer to the vector
            if (itNewW != bfdNew.data.end())
                vecN = &itNewW->second;
            
            // Build TGraphErrors if old has data
            TGraphErrors* grO=nullptr;
            if (vecO && !vecO->empty()) {
                grO = new TGraphErrors((int)vecO->size());
                for (int i=0; i<(int)vecO->size(); i++){
                    grO->SetPoint(i, (*vecO)[i].x, (*vecO)[i].y);
                    grO->SetPointError(i, 0., (*vecO)[i].e);
                }
                // Assign color, marker, etc.
                grO->SetMarkerStyle(20);
                grO->SetMarkerColor(cols[idx%8]);
                grO->SetLineColor(cols[idx%8]);
                grO->SetNameTitle(Form("oldW%d",w), Form("Old(W=%d)",w));
                mgOld->Add(grO, "P");
            }
            
            // Build TGraphErrors if new has data
            TGraphErrors* grN=nullptr;
            if (vecN && !vecN->empty()) {
                grN = new TGraphErrors((int)vecN->size());
                for (int i=0; i<(int)vecN->size(); i++){
                    grN->SetPoint(i, (*vecN)[i].x, (*vecN)[i].y);
                    grN->SetPointError(i, 0., (*vecN)[i].e);
                }
                // Assign color, marker, etc.
                grN->SetMarkerStyle(21);
                grN->SetMarkerColor(cols[idx%8]);
                grN->SetLineColor(cols[idx%8]);
                grN->SetNameTitle(Form("newW%d",w), Form("New(W=%d)",w));
                mgNew->Add(grN, "P");
            }
            
            // Build ratio TGraph if *both* old & new exist
            if (grO && grN && vecO && vecN)
            {
                // Easiest way => build a map oldVals[ x ] => y,
                // then for each new data point at the same x, ratio = new / old
                std::map<double,double> oldVals;
                for (auto &pt : *vecO)
                    oldVals[ pt.x ] = pt.y;
                
                std::vector<double> rx, ry;
                for (auto &ptN : *vecN) {
                    double dx = ptN.x;
                    if (oldVals.count(dx) > 0 && oldVals[dx] > 0) {
                        rx.push_back(dx);
                        ry.push_back( ptN.y / oldVals[dx] );
                    }
                }
                if (!rx.empty()) {
                    TGraph* grR = new TGraph((int)rx.size());
                    for (int i=0; i<(int)rx.size(); i++){
                        grR->SetPoint(i, rx[i], ry[i]);
                    }
                    grR->SetMarkerStyle(25);
                    grR->SetMarkerColor(cols[idx%8]);
                    grR->SetLineColor(cols[idx%8]);
                    grR->SetNameTitle(Form("ratioW%d", w), Form("Ratio(W=%d)", w));
                    mgRatio->Add(grR, "P");
                }
            }
            
            idx++;
        } // end wSet
        
        // The “bestW” for OLD is bfdOld.bestW; for NEW is bfdNew.bestW
        int bestW_old = bfdOld.bestW;
        int bestW_new = bfdNew.bestW;
        
        // Build a 2×3 canvas => same layout as your original
        TCanvas cCmp(Form("compareAll_s%d_ib%d", s,b), "", 2000,2400);
        cCmp.Divide(2,3,0.001, 0.001);
        
        // 1) top-left => mgOld
        cCmp.cd(1);
        {
            int nGo = (mgOld->GetListOfGraphs() ? mgOld->GetListOfGraphs()->GetSize() : 0);
            if (nGo>0) {
                mgOld->Draw("A");
                mgOld->GetXaxis()->SetLimits(-2000,2000);
                mgOld->SetMinimum(0);
                mgOld->SetMaximum(16000);
                gPad->SetGrid();
            } else {
                TLatex la; la.DrawLatexNDC(0.3,0.5,"No OLD data!");
            }
        }
        
        // 2) top-right => mgNew
        cCmp.cd(2);
        {
            int nGn = (mgNew->GetListOfGraphs() ? mgNew->GetListOfGraphs()->GetSize() : 0);
            if (nGn>0) {
                mgNew->Draw("A");
                mgNew->GetXaxis()->SetLimits(-2000,2000);
                mgNew->SetMinimum(0);
                mgNew->SetMaximum(16000);
                gPad->SetGrid();
            } else {
                TLatex la; la.DrawLatexNDC(0.3,0.5,"No NEW data!");
            }
        }
        
        // 3) middle-left => ratio
        cCmp.cd(3);
        gPad->SetGrid();
        {
            int nr = (mgRatio->GetListOfGraphs() ? mgRatio->GetListOfGraphs()->GetSize() : 0);
            if (nr>0) {
                mgRatio->Draw("A");
            } else {
                TLatex la; la.DrawLatexNDC(0.3,0.5,"No ratio data!");
            }
        }
        
        // 4) middle-right => legend
        cCmp.cd(4);
        TLegend* leg = new TLegend(0.1,0.1,0.9,0.9);
        leg->SetNColumns(3);
        leg->SetFillColorAlpha(kWhite,0.8);
        leg->SetBorderSize(1);
        leg->SetTextSize(0.03);
        
        auto addLegend = [&](TMultiGraph* mg){
            if(!mg->GetListOfGraphs()) return;
            for (int i=0; i < mg->GetListOfGraphs()->GetSize(); i++){
                auto gObj = mg->GetListOfGraphs()->At(i);
                if (auto tg = dynamic_cast<TGraph*>(gObj)) {
                    leg->AddEntry(tg, tg->GetTitle(), "lp");
                }
            }
        };
        addLegend(mgOld);
        addLegend(mgNew);
        addLegend(mgRatio);
        leg->Draw();
        
        // ----------------------------------------------------------------------
        // bottom-left => old bestW fit
        // ----------------------------------------------------------------------
        double kOld = 0.;
        double oldA = 0.;
        cCmp.cd(5);
        {
            // bfdOld.data[ bestW_old ] => vector<XYerr>
            auto itData = bfdOld.data.find(bestW_old);
            if (bestW_old<1 || itData==bfdOld.data.end() || itData->second.empty()) {
                TLatex la;
                la.DrawLatexNDC(0.3,0.5,"No OLD bestW data");
            }
            else
            {
                const auto &arr = itData->second;
                TGraphErrors* gE = new TGraphErrors((int)arr.size());
                
                double minY=1e9, maxY=-1e9;
                for (int i=0; i<(int)arr.size(); i++){
                    gE->SetPoint(i, arr[i].x, arr[i].y);
                    gE->SetPointError(i, 0., arr[i].e);
                    if (arr[i].y<minY) minY=arr[i].y;
                    if (arr[i].y>maxY) maxY=arr[i].y;
                }
                double yLo= (minY>0 ? 0.8*minY : 0.);
                double yHi= 1.2*maxY;
                
                gPad->DrawFrame(-2000,yLo,2000,yHi,
                                Form("OLD bestW=%d;#Delta(mV);AvgADC", bestW_old));
                
                gE->SetMarkerStyle(20);
                gE->SetMarkerColor(kBlue+1);
                gE->Draw("P SAME");
                
                // We already have amplitude/k in bfdOld, so we do NOT strictly need to re-fit.
                // But if you want to re-draw the exponential with the stored amplitude/k:
                TF1* fExp = new TF1("fExp_old","[0]*exp([1]*x)", -1500,1500);
                fExp->SetLineColor(kRed);
                fExp->SetParameter(0, bfdOld.amplitude);
                fExp->SetParameter(1, bfdOld.kParam);
                fExp->Draw("SAME");
                
                // For nice prints, store them in local
                oldA  = bfdOld.amplitude;
                kOld  = bfdOld.kParam;
                double eA = bfdOld.amplitudeErr;
                double eK = bfdOld.kParamErr;
                
                // Just for consistency, define a local "resolution(k)"
                double resolutionK=0.;
                if (std::fabs(kOld)>1e-12) {
                    resolutionK= (eK/std::fabs(kOld))*100.;
                }
                
                // Draw latex on canvas
                {
                    TLatex lat;
                    lat.SetNDC(true);
                    lat.SetTextSize(0.04);
                    double yT=0.78;
                    lat.DrawLatex(0.15,yT, Form("A=%.2f #pm %.2f", oldA,eA));
                    yT -= 0.06;
                    lat.DrawLatex(0.15,yT, Form("k=%.4g #pm %.4g", kOld,eK));
                    yT -= 0.06;
                    lat.DrawLatex(0.15,yT, Form("Res(k)= %.2f%%", resolutionK));
                }
            }
        }
        
        // ----------------------------------------------------------------------
        // bottom-right => new bestW fit
        // ----------------------------------------------------------------------
        double kNew = 0.;
        double newA = 0.;
        cCmp.cd(6);
        {
            auto itData = bfdNew.data.find(bestW_new);
            if (bestW_new<1 || itData==bfdNew.data.end() || itData->second.empty()) {
                TLatex la;
                la.DrawLatexNDC(0.3,0.5,"No NEW bestW data");
            }
            else
            {
                const auto &arr = itData->second;
                TGraphErrors* gE = new TGraphErrors((int)arr.size());
                
                double minY=1e9, maxY=-1e9;
                for (int i=0; i<(int)arr.size(); i++){
                    gE->SetPoint(i, arr[i].x, arr[i].y);
                    gE->SetPointError(i, 0., arr[i].e);
                    if (arr[i].y<minY) minY=arr[i].y;
                    if (arr[i].y>maxY) maxY=arr[i].y;
                }
                double yLo= (minY>0 ? 0.8*minY : 0.);
                double yHi= 1.2*maxY;
                
                gPad->DrawFrame(-2000,yLo,2000,yHi,
                                Form("NEW bestW=%d;#Delta(mV);AvgADC", bestW_new));
                
                gE->SetMarkerStyle(21);
                gE->SetMarkerColor(kRed+1);
                gE->Draw("P SAME");
                
                TF1* fExp = new TF1("fExp_new","[0]*exp([1]*x)", -1500,1500);
                fExp->SetLineColor(kBlue);
                fExp->SetParameter(0, bfdNew.amplitude);
                fExp->SetParameter(1, bfdNew.kParam);
                fExp->Draw("SAME");
                
                newA = bfdNew.amplitude;
                kNew = bfdNew.kParam;
                double eA = bfdNew.amplitudeErr;
                double eK = bfdNew.kParamErr;
                double resolutionK=0.;
                if (std::fabs(kNew)>1e-12) {
                    resolutionK= (eK/std::fabs(kNew))*100.;
                }
                
                // draw latex
                {
                    TLatex lat;
                    lat.SetNDC(true);
                    lat.SetTextSize(0.04);
                    double yT=0.78;
                    lat.DrawLatex(0.15,yT, Form("A=%.2f #pm %.2f", newA,eA));
                    yT -= 0.06;
                    lat.DrawLatex(0.15,yT, Form("k=%.4g #pm %.4g", kNew,eK));
                    yT -= 0.06;
                    lat.DrawLatex(0.15,yT, Form("Res(k)= %.2f%%", resolutionK));
                }
            }
        }
        
        // 5) Store final (k) + amplitude in your existing arrays
        kRatios.push_back({ s,b, kOld, kNew });
        amplitudeRatios.push_back({ s,b, oldA, newA });
        
        // Save final multiW overlay
        std::string outName = Form("%s/compareAll_s%d_ib%d.png", multiW.c_str(), s,b);
        cCmp.SaveAs(outName.c_str());
    } // end loop over oldFitMap
    
    std::cout << "[INFO] Created " << validIB
              << " IB multiW overlay plots in " << multiW << "\n\n";
    
    // Additional summary: like your final block that prints min..max.. etc.
    {
        double sumKold=0., sumKnew=0.;
        double minKold=1e9, maxKold=-1e9;
        double minKnew=1e9, maxKnew=-1e9;
        int count=0;
        
        for (auto &kk : kRatios) {
            // only consider positive k
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

void produceCompare32vsBest(
    const std::map< std::pair<int,int>, BoardFitData > &oldFitMap,
    const std::map< std::pair<int,int>, BoardFitData > &newFitMap,
    const std::string &summaryOutputDir
)
{
    std::cout << "\n[DEBUG] produceCompare32vsBest() => building hist: old(32,0) vs old(best), new(32,0) vs new(best)\n";

    //----------------------------------------------------------------------
    // 1) Create the four histograms:
    //----------------------------------------------------------------------
    TH1F* histOld32   = new TH1F("histOld32",   "OLD=width32, offset=0;ADC;Count",128,0,16384);
    histOld32->SetLineColor(kBlue);
    histOld32->SetLineWidth(3);

    TH1F* histOldBest = new TH1F("histOldBest","OLD bestW offset=0;ADC;Count",128,0,16384);
    histOldBest->SetLineColor(kMagenta+2);
    histOldBest->SetLineWidth(3);

    TH1F* histNew32   = new TH1F("histNew32","NEW=width32, offset=0;ADC;Count",128,0,16384);
    histNew32->SetLineColor(kBlue);
    histNew32->SetLineWidth(3);

    TH1F* histNewBest = new TH1F("histNewBest","NEW bestW offset=0;ADC;Count",128,0,16384);
    histNewBest->SetLineColor(kRed);
    histNewBest->SetLineWidth(3);

    // We'll also track newBestContributors => bin -> vector of (s,b),
    // so we can print out the boards that cause spikes in [2500..3500].
    std::map<int, std::vector< std::pair<int,int> >> newBestContributors;

    //----------------------------------------------------------------------
    // 2) Fill histograms for OLD sample => (width=32) and (bestW)
    //----------------------------------------------------------------------
    for (auto &pr : oldFitMap)
    {
        const BoardFitData &bfd = pr.second;

        // (A) bestW
        int wBest = bfd.bestW;
        if (wBest > 0) {
            // bfd.data[wBest] => vector<XYerr>
            auto itVec = bfd.data.find(wBest);
            if (itVec != bfd.data.end()) {
                for (auto &xy : itVec->second) {
                    if (std::fabs(xy.x)<1e-9 && xy.y>0.) {
                        histOldBest->Fill(xy.y);
                    }
                }
            }
        }

        // (B) width=32
        auto it32 = bfd.data.find(32);
        if (it32 != bfd.data.end()) {
            for (auto &xy : it32->second) {
                if (std::fabs(xy.x)<1e-9 && xy.y>0.) {
                    histOld32->Fill(xy.y);
                }
            }
        }
    }

    //----------------------------------------------------------------------
    // 3) Fill histograms for NEW sample => (width=32) and (bestW)
    //----------------------------------------------------------------------
    for (auto &pr : newFitMap)
    {
        // (A) bestW
        const BoardFitData &bfd = pr.second;
        int s = bfd.sector;
        int b = bfd.ib;
        int wBest = bfd.bestW;

        if (wBest > 0) {
            auto itVec = bfd.data.find(wBest);
            if (itVec != bfd.data.end()) {
                for (auto &xy : itVec->second) {
                    if (std::fabs(xy.x)<1e-9 && xy.y>0.) {
                        double adcVal = xy.y;
                        histNewBest->Fill(adcVal);

                        // track => which boards contributed to this bin
                        int binIndex = histNewBest->FindBin(adcVal);
                        newBestContributors[binIndex].push_back({ s,b });
                    }
                }
            }
        }

        // (B) width=32
        auto it32 = bfd.data.find(32);
        if (it32 != bfd.data.end()) {
            for (auto &xy : it32->second) {
                if (std::fabs(xy.x)<1e-9 && xy.y>0.) {
                    histNew32->Fill(xy.y);
                }
            }
        }
    }

    //----------------------------------------------------------------------
    // 4) Print out boards that cause “spikes” in NEW bestW, ADC in [2500..3500]
    //----------------------------------------------------------------------
    std::cout << "\n[DEBUG] => Listing boards that cause spikes in NEW bestW, ADC in [2500..3500]...\n";
    int firstBin = histNewBest->FindBin(2500.);
    int lastBin  = histNewBest->FindBin(3500.);
    for(int ibin = firstBin; ibin <= lastBin; ibin++){
        int count = (int) histNewBest->GetBinContent(ibin);
        if (count<1) continue;

        double center = histNewBest->GetBinCenter(ibin);
        std::cout << "  Bin #" << ibin << " => center=" << center
                  << ", count=" << count << ":\n";

        auto &vec = newBestContributors[ibin];
        for (auto &sb : vec) {
            int sec = sb.first;
            int ibd = sb.second;
            std::cout << "     sector=" << sec << ", IB=" << ibd << "\n";
        }
    }

    //----------------------------------------------------------------------
    // 5) A small helper => binStats(...) returns #below, #inRange, #above, #total
    //----------------------------------------------------------------------
    auto binStats = [&](TH1F *hist){
        int below=0, inRange=0, above=0;
        int nbins= hist->GetNbinsX();
        for(int i=1; i<=nbins; i++){
            double xC   = hist->GetBinCenter(i);
            int    c    = (int)hist->GetBinContent(i);
            if(xC < 2500.) below   += c;
            else if(xC > 3500.) above += c;
            else inRange += c;
        }
        int total = below + inRange + above;
        return std::make_tuple(below, inRange, above, total);
    };

    //----------------------------------------------------------------------
    // 6) Build the final 2×1 TCanvas => "compareADC_32_vs_best.png"
    //    with top sub‐pads showing the hist + overlaid lines,
    //    and bottom sub‐pads for ratio hist.
    //----------------------------------------------------------------------
    TCanvas cF("cF","compareADC_32_vs_best(OLD vs NEW)",1800,700);
    cF.Divide(2,1);

    //--------------------------------------------------------------------
    // LEFT => old(32 vs best)
    //--------------------------------------------------------------------
    cF.cd(1);
    TPad* subPad1= (TPad*)gPad;
    subPad1->Divide(1,2,0.001,0.001);

    // top => old
    subPad1->cd(1);
    double mRef  = histOld32->GetMaximum();
    double mBest = histOldBest->GetMaximum();
    double yMax  = std::max(mRef,mBest) * 1.2;

    histOld32->SetMaximum(yMax);
    histOld32->SetTitle("OLD Data: width=32 vs. best (#Delta=0)");
    histOld32->GetXaxis()->SetTitle("Average ADC");
    histOld32->GetYaxis()->SetTitle("Count of IBs");
    histOld32->Draw("HIST");
    histOldBest->Draw("HIST SAME");

    // green lines at [2500..3500]
    TLine lA(2500,0,2500,yMax);
    lA.SetLineColor(kGreen+2); lA.SetLineStyle(2); lA.SetLineWidth(2);
    lA.Draw("same");
    TLine lB(3500,0,3500,yMax);
    lB.SetLineColor(kGreen+2); lB.SetLineStyle(2); lB.SetLineWidth(2);
    lB.Draw("same");

    TLegend legO(0.42,0.70,0.77,0.88);
    legO.AddEntry(histOld32,   "OLD Sample (TP=32)",   "l");
    legO.AddEntry(histOldBest, "OLD Sample (best TP)", "l");
    legO.AddEntry(&lA,         "Range: [2500..3500]",  "l");
    legO.Draw();

    // print tallies
    {
        auto [b32,i32,a32,t32]        = binStats(histOld32);
        auto [bBest,iBest,aBest,tBest]= binStats(histOldBest);

        std::cout << "\n[DEBUG] OLD(32,0) => total IB="<< t32
                  << " => below="<< b32 <<", inRange="<< i32 <<", above="<< a32 <<"\n";
        if(t32>0){
            std::cout << "         => % below= "<< (100.*b32/t32)
                      <<", % inRange= "<< (100.*i32/t32)
                      <<", % above= "<< (100.*a32/t32) <<"\n";
        }
        std::cout << "[DEBUG] OLD(best,0) => total IB="<< tBest
                  << " => below="<< bBest <<", inRange="<< iBest <<", above="<< aBest <<"\n";
        if(tBest>0){
            std::cout << "         => % below= "<< (100.*bBest/tBest)
                      <<", % inRange= "<< (100.*iBest/tBest)
                      <<", % above= "<< (100.*aBest/tBest) <<"\n";
        }

        // small text blocks
        TLatex latL, latR;
        latL.SetNDC(true); latR.SetNDC(true);
        latL.SetTextSize(0.04);
        latR.SetTextSize(0.04);

        double yLeft = 0.6;
        latL.SetTextColor(kBlack);
        latL.DrawLatex(0.35,yLeft, Form("OLD(32) total IB=%d", t32));
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
        latR.DrawLatex(0.63,yRight, Form("OLD(best) total IB=%d", tBest));
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

    // bottom => ratio
    subPad1->cd(2);
    TH1F* hRatioO= (TH1F*) histOldBest->Clone("ratioOld");
    hRatioO->Divide(histOld32);
    hRatioO->SetTitle("");
    hRatioO->GetXaxis()->SetTitle("ADC");
    hRatioO->GetYaxis()->SetTitle("OLD(best)/OLD(32)");
    hRatioO->GetYaxis()->SetNdivisions(505);

    double rmO=0.;
    for(int i=1; i<= hRatioO->GetNbinsX(); i++){
        double val= hRatioO->GetBinContent(i);
        double err= hRatioO->GetBinError(i);
        if(val+err> rmO) rmO= val+err;
    }
    hRatioO->SetMaximum(std::max(3.0,1.2*rmO));
    hRatioO->SetMinimum(0.);
    hRatioO->Draw("E0");

    //--------------------------------------------------------------------
    // RIGHT => new(32 vs best)
    //--------------------------------------------------------------------
    cF.cd(2);
    TPad* subPad2= (TPad*)gPad;
    subPad2->Divide(1,2,0.001,0.001);

    // top => new
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
    LN1.SetLineColor(kGreen+2); LN1.SetLineStyle(2); LN1.SetLineWidth(2);
    LN1.Draw("same");
    TLine LN2(3500,0,3500,yM);
    LN2.SetLineColor(kGreen+2); LN2.SetLineStyle(2); LN2.SetLineWidth(2);
    LN2.Draw("same");

    TLegend legN(0.42,0.70,0.77,0.88);
    legN.AddEntry(histNew32,   "NEW Sample (TP=32)",   "l");
    legN.AddEntry(histNewBest, "NEW Sample (best TP)", "l");
    legN.AddEntry(&LN1,        "Range: [2500..3500]",  "l");
    legN.Draw();

    // print tallies
    {
        auto [b32,i32,a32,t32]        = binStats(histNew32);
        auto [bBest,iBest,aBest,tBest]= binStats(histNewBest);

        std::cout << "\n[DEBUG] NEW(32,0) => total IB="<< t32
                  << " => below="<< b32 <<", inRange="<< i32 <<", above="<< a32 <<"\n";
        if(t32>0){
            std::cout << "         => % below= "<< (100.*b32/t32)
                      <<", % inRange= "<< (100.*i32/t32)
                      <<", % above= "<< (100.*a32/t32) <<"\n";
        }
        std::cout << "[DEBUG] NEW(best,0) => total IB="<< tBest
                  << " => below="<< bBest <<", inRange="<< iBest <<", above="<< aBest <<"\n";
        if(tBest>0){
            std::cout << "         => % below= "<< (100.*bBest/tBest)
                      <<", % inRange= "<< (100.*iBest/tBest)
                      <<", % above= "<< (100.*aBest/tBest) <<"\n";
        }

        TLatex latL, latR;
        latL.SetNDC(true); latR.SetNDC(true);
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

    // bottom => ratio
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
        if(val+err> rmN) rmN= val+err;
    }
    hRatioN->SetMaximum(std::max(3.0,1.2*rmN));
    hRatioN->SetMinimum(0.);
    hRatioN->Draw("E0");

    //----------------------------------------------------------------------
    // 7) Save final figure
    //----------------------------------------------------------------------
    std::string outTwo = summaryOutputDir + "/compareADC_32_vs_best.png";
    cF.SaveAs(outTwo.c_str());
    std::cout << "[INFO] => Created " << outTwo << " via produceCompare32vsBest()\n";
}


void produceKAndAmplitudePlots(
    const std::map< std::pair<int,int>, BoardFitData > &oldFitMap,
    const std::map< std::pair<int,int>, BoardFitData > &newFitMap,
    const std::string &summary
)
{
    // 1) Build kRatios + amplitudeRatios by reading (kParam, amplitude).
    std::vector<Kval>   kRatios;
    std::vector<AmpVal> amplitudeRatios;

    for (auto &pr : oldFitMap) {
        int s = pr.first.first;
        int b = pr.first.second;

        // Check if newFitMap has this same (sector, ib)
        auto itNew = newFitMap.find({s,b});
        if (itNew == newFitMap.end()) {
            // not found => skip
            continue;
        }

        const BoardFitData &bfdOld = pr.second;
        const BoardFitData &bfdNew = itNew->second;

        double kO   = bfdOld.kParam;
        double kN   = bfdNew.kParam;
        double ampO = bfdOld.amplitude;
        double ampN = bfdNew.amplitude;

        // Store them
        kRatios.push_back(  {s,b, kO,   kN} );
        amplitudeRatios.push_back( {s,b, ampO, ampN} );
    }

    //----------------------------------------------------------------------
    // Step7 => produce k ratio => TGraph => k(new)/k(old)
    //----------------------------------------------------------------------
    std::cout << "[DEBUG] => produceKAndAmplitudePlots() => Step7 => building k ratio plot...\n";
    {
        TGraph* grK = new TGraph();
        grK->SetMarkerStyle(20);
        grK->SetMarkerColor(kBlue+1);
        grK->SetLineColor(kBlue+1);

        int iPt=0;
        for (auto &kk : kRatios) {
            double ratio=0.;
            if (kk.kOld>0.) {
                ratio= kk.kNew / kk.kOld;
            }
            grK->SetPoint(iPt, iPt+1, ratio);
            iPt++;
        }

        int nPts= iPt;
        if (nPts < 1) {
            std::cout << "[INFO] => No boards => no k ratio plot.\n";
            delete grK;
        } else {
            TCanvas cK("cK","k ratio new/old",900,600);
            double xMax= nPts + 1.;
            double rm=1.;
            // find the max ratio for setting y-axis
            for(int i=0; i<nPts; i++){
                double xx,yy;
                grK->GetPoint(i, xx, yy);
                if (yy>rm) rm= yy;
            }

            TH1F* frm= cK.DrawFrame(0., 0., xMax, 1.2*rm);
            frm->SetTitle("k ratio(new/old) for bestW fits;Index; ratio");
            grK->Draw("P SAME");

            std::string outK = summary + "/kRatio_bestW.png";
            cK.SaveAs(outK.c_str());
            std::cout << "[INFO] => Created " << outK << "\n";

            delete grK;
        }
    }

    //----------------------------------------------------------------------
    // Step8 => side-by-side 1D distributions => hKold, hKnew => gauss fits
    //----------------------------------------------------------------------
    // (kept EXACTLY as your snippet requires)
    //----------------------------------------------------------------------
    std::cout << "[DEBUG] => produceKAndAmplitudePlots() => Step8 => side-by-side 1D distributions for k...\n";
    {
        double maxKold=0., maxKnew=0.;
        int countOld=0, countNew=0;

        // find maxK
        for (auto &kk : kRatios) {
            if (kk.kOld > maxKold) maxKold= kk.kOld;
            if (kk.kNew > maxKnew) maxKnew= kk.kNew;
        }

        double maxK= std::max(maxKold, maxKnew);
        if (maxK<1e-7) maxK= 1e-5;  // avoid 0 range
        double xMax= 1.2 * maxK;

        // create hist
        TH1F* hKold= new TH1F("hKold","Old k distribution; k; # IBs",50,0., xMax);
        hKold->SetLineColor(kBlue+2);
        hKold->SetLineWidth(3);

        TH1F* hKnew= new TH1F("hKnew","New k distribution; k; # IBs",50,0., xMax);
        hKnew->SetLineColor(kRed+1);
        hKnew->SetLineWidth(3);

        // fill
        for (auto &kk : kRatios) {
            if (kk.kOld>0.) { hKold->Fill(kk.kOld); countOld++; }
            if (kk.kNew>0.) { hKnew->Fill(kk.kNew); countNew++; }
        }

        // debug prints for suspicious
        for (auto &kk : kRatios) {
            // define "suspicious" as k < 1e-6 or k > 1e-3
            if (kk.kOld>0. && (kk.kOld<1e-6 || kk.kOld>1e-3)) {
                std::cout << "[DEBUG] *** Outlier in OLD *** => s="<<kk.sector
                          <<", ib="<<kk.ib<<", kOld="<<kk.kOld<<"\n";
            }
            if (kk.kNew>0. && (kk.kNew<1e-6 || kk.kNew>1e-3)) {
                std::cout << "[DEBUG] *** Outlier in NEW *** => s="<<kk.sector
                          <<", ib="<<kk.ib<<", kNew="<<kk.kNew<<"\n";
            }
        }

        TCanvas cKD("cKD","k distributions",1200,600);
        cKD.Divide(2,1);

        // left => old
        cKD.cd(1);
        hKold->Draw("HIST");
        TF1* fGold= new TF1("fGold","gaus",0., xMax);
        fGold->SetLineColor(kGreen+3);
        fGold->SetLineWidth(2);
        fGold->SetNpx(2000);
        hKold->Fit(fGold,"Q");
        fGold->Draw("SAME");

        double meanOld   = fGold->GetParameter(1);
        double eMeanOld  = fGold->GetParError(1);
        double sigOld    = fGold->GetParameter(2);
        double eSigOld   = fGold->GetParError(2);

        double resOld=0.;
        if (std::fabs(meanOld)>1e-12) {
            resOld= sigOld / std::fabs(meanOld);
        }

        {
            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextSize(0.04);
            double yT=0.80;
            lat.DrawLatex(0.50,yT,Form("N_{IB} = %d",countOld));          yT-=0.06;
            lat.DrawLatex(0.50,yT,Form("#mu=%.3g #pm %.3g",meanOld,eMeanOld));  yT-=0.06;
            lat.DrawLatex(0.50,yT,Form("#sigma=%.3g #pm %.3g",sigOld,eSigOld)); yT-=0.06;
            lat.DrawLatex(0.50,yT,Form("Res=%.1f%%",resOld*100.));
        }

        // right => new
        cKD.cd(2);
        hKnew->Draw("HIST");
        TF1* fGnew= new TF1("fGnew","gaus",0., xMax);
        fGnew->SetLineColor(kGreen+3);
        fGnew->SetLineWidth(2);
        fGnew->SetNpx(2000);
        hKnew->Fit(fGnew,"Q");
        fGnew->Draw("SAME");

        double meanNew   = fGnew->GetParameter(1);
        double eMeanNew  = fGnew->GetParError(1);
        double sigNew    = fGnew->GetParameter(2);
        double eSigNew   = fGnew->GetParError(2);

        double resNew=0.;
        if (std::fabs(meanNew)>1e-12) {
            resNew= sigNew / std::fabs(meanNew);
        }

        {
            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextSize(0.04);
            double yT=0.80;
            lat.DrawLatex(0.50,yT,Form("N_{IB} = %d",countNew));          yT-=0.06;
            lat.DrawLatex(0.50,yT,Form("#mu=%.3g #pm %.3g",meanNew,eMeanNew));  yT-=0.06;
            lat.DrawLatex(0.50,yT,Form("#sigma=%.3g #pm %.3g",sigNew,eSigNew)); yT-=0.06;
            lat.DrawLatex(0.50,yT,Form("Res=%.1f%%",resNew*100.));
        }

        std::string outKdist = summary + "/kValue_distributions.png";
        cKD.SaveAs(outKdist.c_str());
        std::cout << "[INFO] => Created " << outKdist << "\n";

        // cleanup
        delete fGold;
        delete fGnew;
        delete hKold;
        delete hKnew;
    }

    //----------------------------------------------------------------------
    // Step10 => side-by-side amplitude distributions (Landau) + ratio
    //----------------------------------------------------------------------
    std::cout << "[DEBUG] => produceKAndAmplitudePlots() => amplitude distributions + ratio...\n";
    {
        // 1) side-by-side amplitude distributions => Landau fits
        double maxAold=0., maxAnew=0.;
        int countAold=0, countAnew=0;
        for (auto &av : amplitudeRatios) {
            if (av.ampOld> maxAold) maxAold= av.ampOld;
            if (av.ampNew> maxAnew) maxAnew= av.ampNew;
        }

        double maxA= std::max(maxAold, maxAnew);
        if (maxA<1e-7) maxA= 1e-5;  // avoid zero range
        double xMax= 1.2 * maxA;

        TH1F* hAold= new TH1F("hAold","Old amplitude distribution (Landau); A; # IBs",
                              50, 0., xMax);
        hAold->SetLineColor(kBlue+2);
        hAold->SetLineWidth(3);

        TH1F* hAnew= new TH1F("hAnew","New amplitude distribution (Landau); A; # IBs",
                              50, 0., xMax);
        hAnew->SetLineColor(kRed+1);
        hAnew->SetLineWidth(3);

        // fill
        for (auto &av : amplitudeRatios) {
            if (av.ampOld>0.) { hAold->Fill(av.ampOld); countAold++; }
            if (av.ampNew>0.) { hAnew->Fill(av.ampNew); countAnew++; }
        }

        TCanvas cAD("cAD","Amplitude distributions (Landau)",1200,600);
        cAD.Divide(2,1);

        // left => old amplitude
        cAD.cd(1);
        hAold->Draw("HIST");

        TF1* fAold= new TF1("fAold","landau", 0., xMax);
        fAold->SetLineColor(kGreen+3);
        fAold->SetLineWidth(2);
        fAold->SetNpx(2000);

        hAold->Fit(fAold,"Q");
        fAold->Draw("SAME");

        double mpvOld    = fAold->GetParameter(1);
        double eMPVOld   = fAold->GetParError(1);
        double widthOld  = fAold->GetParameter(2);
        double eWidthOld = fAold->GetParError(2);

        // define "resolution" => width / MPV
        double resAold=0.;
        if (std::fabs(mpvOld)>1e-12) {
            resAold= widthOld / std::fabs(mpvOld);
        }

        {
            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextSize(0.04);
            double yT=0.80;
            lat.DrawLatex(0.50,yT,Form("N_{IB} = %d",countAold));           yT-=0.06;
            lat.DrawLatex(0.50,yT,Form("MPV=%.3g #pm %.3g", mpvOld,eMPVOld));   yT-=0.06;
            lat.DrawLatex(0.50,yT,Form("Width=%.3g #pm %.3g", widthOld,eWidthOld)); yT-=0.06;
            lat.DrawLatex(0.50,yT,Form("Res=%.1f%%",resAold*100.));
        }

        // right => new amplitude
        cAD.cd(2);
        hAnew->Draw("HIST");

        TF1* fAnew= new TF1("fAnew","landau", 0., xMax);
        fAnew->SetLineColor(kGreen+3);
        fAnew->SetLineWidth(2);
        fAnew->SetNpx(2000);

        hAnew->Fit(fAnew,"Q");
        fAnew->Draw("SAME");

        double mpvNew    = fAnew->GetParameter(1);
        double eMPVNew   = fAnew->GetParError(1);
        double widthNew  = fAnew->GetParameter(2);
        double eWidthNew = fAnew->GetParError(2);

        double resAnew=0.;
        if (std::fabs(mpvNew)>1e-12) {
            resAnew= widthNew / std::fabs(mpvNew);
        }

        {
            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextSize(0.04);
            double yT=0.80;
            lat.DrawLatex(0.50,yT,Form("N_{IB} = %d",countAnew));           yT-=0.06;
            lat.DrawLatex(0.50,yT,Form("MPV=%.3g #pm %.3g", mpvNew,eMPVNew));    yT-=0.06;
            lat.DrawLatex(0.50,yT,Form("Width=%.3g #pm %.3g", widthNew,eWidthNew)); yT-=0.06;
            lat.DrawLatex(0.50,yT,Form("Res=%.1f%%",resAnew*100.));
        }

        std::string outAmpDist= summary + "/amplitude_distributions_landau.png";
        cAD.SaveAs(outAmpDist.c_str());
        std::cout << "[INFO] => Created " << outAmpDist << "\n";

        // cleanup
        delete fAold;
        delete fAnew;
        delete hAold;
        delete hAnew;
    }

    // amplitude ratio => A(new)/A(old)
    std::cout << "[DEBUG] => produceKAndAmplitudePlots() => amplitude ratio plot...\n";
    {
        TGraph* grA = new TGraph();
        grA->SetMarkerStyle(20);
        grA->SetMarkerColor(kBlue+1);
        grA->SetLineColor(kBlue+1);

        int iPt=0;
        for (auto &av : amplitudeRatios) {
            double ratioA=0.;
            if (av.ampOld>0.) {
                ratioA= av.ampNew / av.ampOld;
            }
            grA->SetPoint(iPt, iPt+1, ratioA);
            iPt++;
        }

        int nPts= iPt;
        if (nPts<1) {
            std::cout << "[INFO] => No boards => no amplitude ratio plot.\n";
            delete grA;
        } else {
            TCanvas cA("cA","Amplitude ratio new/old",900,600);

            double xMax= nPts + 1.;
            double rm=1.;
            for(int i=0; i<nPts; i++){
                double xx,yy;
                grA->GetPoint(i, xx, yy);
                if (yy>rm) rm= yy;
            }

            TH1F* frmA= cA.DrawFrame(0., 0., xMax, 1.2*rm);
            frmA->SetTitle("Amplitude ratio (new/old) for bestW fits;Index; ratio");
            grA->Draw("P SAME");

            std::string outA= summary + "/amplitudeRatio_bestW.png";
            cA.SaveAs(outA.c_str());
            std::cout << "[INFO] => Created " << outA << "\n";

            delete grA;
        }
    }

    // done
    std::cout<<"\n[INFO] => produceKAndAmplitudePlots() => done.\n\n";
}



void produceColorCodedMapsAndExtra(
    const std::map< std::pair<int,int>, BoardFitData > &oldFitMap,
    const std::map< std::pair<int,int>, BoardFitData > &newFitMap,
    const std::map< std::pair<int,int>, bool > &sectorIBmap,
    const std::string &summary
)
{
    //----------------------------------------------------------------------
    // 1) Helper to retrieve the “offset=0” ADC from W=32
    //    (If missing or not found => return -1.)
    //----------------------------------------------------------------------
    auto get32offset0 = [&](bool oldNotNew, int sector, int ib) -> double {
        // For old => oldFitMap, for new => newFitMap
        const auto &mapRef = (oldNotNew ? oldFitMap : newFitMap);

        // If that board doesn’t exist in the map => no data
        auto itB = mapRef.find({sector, ib});
        if (itB == mapRef.end()) return -1.0;

        // If we do have that board => check if data[32] is present
        const BoardFitData &bfd = itB->second;
        auto itW = bfd.data.find(32);
        if (itW == bfd.data.end()) return -1.0;

        // Look for x=0 => offset=0 => pick y>0
        for (auto &xy : itW->second) {
            if (std::fabs(xy.x)<1e-9 && xy.y>0.) {
                return xy.y;
            }
        }
        return -1.0;
    };

    //----------------------------------------------------------------------
    // 2) Helper to retrieve the “offset=0” ADC from bestW
    //----------------------------------------------------------------------
    auto getBestOffset0 = [&](bool oldNotNew, int sector, int ib) -> double {
        // For old => oldFitMap, for new => newFitMap
        const auto &mapRef = (oldNotNew ? oldFitMap : newFitMap);

        auto itB = mapRef.find({sector, ib});
        if (itB == mapRef.end()) return -1.0;

        const BoardFitData &bfd = itB->second;
        int bestW = bfd.bestW;
        if (bestW < 1) return -1.0;  // invalid bestW

        // Check data[ bestW ]
        auto itW = bfd.data.find(bestW);
        if (itW == bfd.data.end()) return -1.0;

        for (auto &xy : itW->second) {
            if (std::fabs(xy.x)<1e-9 && xy.y>0.) {
                return xy.y;
            }
        }
        return -1.0;
    };

    //----------------------------------------------------------------------
    // We'll keep counters for color-coded categories
    //   old_below, old_inrange, old_above, old_none, ...
    // to print at the end.
    //----------------------------------------------------------------------
    int old_below=0, old_inrange=0, old_above=0, old_none=0;
    int new_below=0, new_inrange=0, new_above=0, new_none=0;

    //----------------------------------------------------------------------
    // A helper to get corners for TBox in your tower-based display
    //----------------------------------------------------------------------
    auto getPadCorners = [&](int s, int b){
        int localSec = (s % 32);
        double phiMin= 8.*localSec;
        double phiMax= phiMin + 8.;

        double etaMin=0., etaMax=0.;
        if (s<32) {
            etaMin= 48. + 8.*b;
            etaMax= etaMin + 8.;
        } else {
            etaMin= 40. - 8.*b;
            etaMax= etaMin + 8.;
        }
        if (etaMin>etaMax) std::swap(etaMin, etaMax);
        return std::make_tuple(phiMin, phiMax, etaMin, etaMax);
    };

    //----------------------------------------------------------------------
    // chooseColor => decides which TBox fill color to use, updates counters
    //----------------------------------------------------------------------
    auto chooseColor = [&](double valADC, bool isOld)-> Color_t {
        if (valADC<0.) {
            if(isOld) old_none++; else new_none++;
            return (Color_t)(kGray+1);
        }
        else if (valADC<2500.) {
            if(isOld) old_below++; else new_below++;
            return kBlue;
        }
        else if (valADC>3500.) {
            if(isOld) old_above++; else new_above++;
            return kRed;
        }
        else {
            if(isOld) old_inrange++; else new_inrange++;
            return (Color_t)(kGreen+2);
        }
    };

    //----------------------------------------------------------------------
    // drawColorMapInPad(...) => draws either W=32 or bestW in the *current* pad
    //----------------------------------------------------------------------
    auto drawColorMapInPad = [&](bool oldNotNew, bool w32OrBest){
        // For each known board => pick ADC => fill TBox
        for (auto &sb : sectorIBmap) {
            int s= sb.first.first;
            int b= sb.first.second;
            // pick either W=32 or bestW
            double valADC= (w32OrBest
                            ? get32offset0(oldNotNew, s, b)
                            : getBestOffset0(oldNotNew, s, b));
            auto [x1,x2,y1,y2] = getPadCorners(s,b);

            TBox *box= new TBox(x1,y1, x2,y2);
            box->SetLineColor(kBlack);
            box->SetLineWidth(1);
            Color_t col = chooseColor(valADC, oldNotNew);
            box->SetFillColor(col);
            box->Draw();

            // label => "s.ib"
            double cx= 0.5*(x1 + x2);
            double cy= 0.5*(y1 + y2);
            TLatex lat;
            lat.SetTextAlign(22);
            lat.SetTextSize(0.017);
            lat.SetTextColor(kBlack);
            lat.DrawLatex(cx, cy, Form("%d.%d", s,b));
        }
    };

    //----------------------------------------------------------------------
    // produceCombinedMap(...) => builds the 2×1 color-coded map for either
    // old or new => top= W=32, bottom= bestW
    //----------------------------------------------------------------------
    auto produceCombinedMap = [&](bool oldNotNew, const std::string &outFile) {
        // Decide on naming
        std::string cName  = (oldNotNew ? "cMapOldCombined" : "cMapNewCombined");
        std::string cTitle = (oldNotNew
                              ? "Old color map (W=32 on top, bestW bottom)"
                              : "New color map (W=32 on top, bestW bottom)");

        TCanvas cC(cName.c_str(), cTitle.c_str(), 3000,2400);
        cC.Divide(1,2, 0.001,0.001);

        //------------------
        // TOP => W=32
        //------------------
        cC.cd(1);
        gPad->SetRightMargin(0.20);
        TH2F* hFtop= new TH2F("hFrameTop",";#phi index;#eta index",
                              256,0,256, 96,0,96);
        hFtop->SetStats(false);
        hFtop->GetXaxis()->SetNdivisions(32,false);
        hFtop->GetYaxis()->SetNdivisions(12,false);
        hFtop->Draw("COL");

        drawColorMapInPad(oldNotNew, true);

        // single legend
        TLegend* legC= new TLegend(0.82, 0.10, 0.98, 0.90);
        legC->SetHeader( (oldNotNew
                          ? "OLD: top=W=32, bottom=BEST"
                          : "NEW: top=W=32, bottom=BEST"), "C");
        legC->SetBorderSize(0);
        legC->SetFillColorAlpha(kWhite,0.6);
        legC->SetTextFont(42);
        legC->SetTextSize(0.035);
        legC->SetMargin(0.15);

        // dummy TBoxes for color legend
        TBox bBlue(0,0,1,1);  bBlue.SetFillColor(kBlue);
        TBox bGreen(0,0,1,1); bGreen.SetFillColor(kGreen+2);
        TBox bRed(0,0,1,1);   bRed.SetFillColor(kRed);
        TBox bGray(0,0,1,1);  bGray.SetFillColor(kGray+1);

        legC->AddEntry(&bBlue,  "ADC < 2500","f");
        legC->AddEntry(&bGreen, "2500..3500","f");
        legC->AddEntry(&bRed,   "ADC > 3500","f");
        legC->AddEntry(&bGray,  "No Data","f");
        legC->Draw();

        //------------------
        // BOTTOM => bestW
        //------------------
        cC.cd(2);
        gPad->SetRightMargin(0.20);
        TH2F* hFbot= new TH2F("hFrameBot",";#phi index;#eta index",
                              256,0,256, 96,0,96);
        hFbot->SetStats(false);
        hFbot->GetXaxis()->SetNdivisions(32,false);
        hFbot->GetYaxis()->SetNdivisions(12,false);
        hFbot->Draw("COL");

        drawColorMapInPad(oldNotNew, false);

        // same legend, no second needed

        // final
        cC.SaveAs(outFile.c_str());
        std::cout << "[INFO] => Created combined color-coded map: " << outFile << "\n";

        // cleanup
        delete hFtop;
        delete hFbot;
    };

    //----------------------------------------------------------------------
    // 2) Actually produce the old/new combined color-coded PNGs
    //----------------------------------------------------------------------
    std::cout << "[DEBUG] => producing colorMap_old_combined.png & colorMap_new_combined.png\n";
    {
        std::string outOld= summary + "/colorMap_old_combined.png";
        produceCombinedMap(true, outOld);

        std::string outNew= summary + "/colorMap_new_combined.png";
        produceCombinedMap(false, outNew);
    }

    //----------------------------------------------------------------------
    // 3) Build the single 2×1 “COLZ” figure for the new sample => W=32 vs bestW
    //----------------------------------------------------------------------
    {
        TCanvas cH("cH","New Sample Heat Maps (COLZ)", 2200,1600);
        cH.Divide(1,2);

        //-------------------------
        // TOP pad => new(W=32)
        //-------------------------
        cH.cd(1);
        gPad->SetRightMargin(0.20);

        TH2F* h2W32= new TH2F("h2New32","NEW Sample (TPW=32, #Delta offset=0);#phi index;#eta index",
                              256,0,256, 96,0,96);
        h2W32->SetStats(false);
        h2W32->GetZaxis()->SetTitle("avg ADC");

        // Fill from get32offset0(false => "new")
        for (auto &sb : sectorIBmap) {
            int s= sb.first.first;
            int b= sb.first.second;
            double val= get32offset0(false, s,b);
            if (val<0) continue;

            // Fill an 8×8 region
            int localSec= (s % 32);
            int phiMin= 8*localSec, phiMax= phiMin+7;
            int etaMin=0, etaMax=0;
            if(s<32) {
                etaMin=48+ 8*b;
                etaMax=etaMin+7;
            } else {
                int base= 40 - 8*b;
                etaMin= base; etaMax= base+7;
            }
            if(etaMin>etaMax) std::swap(etaMin, etaMax);

            for(int xx=phiMin; xx<=phiMax; xx++){
                for(int yy= etaMin; yy<=etaMax; yy++){
                    h2W32->SetBinContent(xx+1,yy+1, val);
                }
            }
        }
        h2W32->Draw("COLZ");

        // Also label each IB => "s.b"
        for (auto &sb : sectorIBmap) {
            int s= sb.first.first;
            int b= sb.first.second;
            double val= get32offset0(false,s,b);
            if (val<0) continue;

            int localSec= (s % 32);
            int phiMin= 8*localSec, phiMax= phiMin+7;
            int etaMin=0, etaMax=0;
            if(s<32) {
                etaMin=48 + 8*b;
                etaMax=etaMin+7;
            } else {
                int base=40 - 8*b;
                etaMin= base; etaMax= base+7;
            }
            if(etaMin>etaMax) std::swap(etaMin,etaMax);

            double cx= 0.5*(phiMin + phiMax +1);
            double cy= 0.5*(etaMin + etaMax +1);

            TLatex lat;
            lat.SetTextAlign(22);
            lat.SetTextSize(0.015);
            lat.SetTextColor(kBlack);
            lat.DrawLatex(cx,cy, Form("%d.%d",s,b));
        }

        //-------------------------
        // BOTTOM pad => new(bestW)
        //-------------------------
        cH.cd(2);
        gPad->SetRightMargin(0.20);

        TH2F* h2Best= new TH2F("h2NewBest","NEW (Best W, offset=0);#phi index;#eta index",
                               256,0,256, 96,0,96);
        h2Best->SetStats(false);
        h2Best->GetZaxis()->SetTitle("avg ADC");

        // fill from getBestOffset0(false => new)
        for (auto &sb : sectorIBmap) {
            int s= sb.first.first;
            int b= sb.first.second;
            double val= getBestOffset0(false,s,b);
            if(val<0) continue;

            int localSec= (s % 32);
            int phiMin= 8*localSec, phiMax= phiMin+7;
            int etaMin=0, etaMax=0;
            if(s<32) {
                etaMin= 48 + 8*b;
                etaMax= etaMin+7;
            } else {
                int base=40 - 8*b;
                etaMin= base; etaMax= base+7;
            }
            if(etaMin>etaMax) std::swap(etaMin,etaMax);

            for(int xx=phiMin; xx<=phiMax; xx++){
                for(int yy=etaMin; yy<=etaMax; yy++){
                    h2Best->SetBinContent(xx+1,yy+1, val);
                }
            }
        }
        h2Best->Draw("COLZ");

        // label each IB
        for (auto &sb : sectorIBmap) {
            int s= sb.first.first;
            int b= sb.first.second;
            double val= getBestOffset0(false,s,b);
            if(val<0) continue;

            int localSec= (s % 32);
            int phiMin= 8*localSec, phiMax= phiMin+7;
            int etaMin=0, etaMax=0;
            if(s<32) {
                etaMin=48+8*b; etaMax=etaMin+7;
            } else {
                int base=40 - 8*b;
                etaMin= base; etaMax= base+7;
            }
            if(etaMin>etaMax) std::swap(etaMin,etaMax);

            double cx= 0.5*(phiMin + phiMax +1);
            double cy= 0.5*(etaMin + etaMax +1);

            TLatex lat;
            lat.SetTextAlign(22);
            lat.SetTextSize(0.015);
            lat.SetTextColor(kBlack);
            lat.DrawLatex(cx,cy, Form("%d.%d",s,b));
        }

        // final
        std::string outC= summary + "/heatMap_new_combined_colz.png";
        cH.SaveAs(outC.c_str());
        std::cout<<"[INFO] => Created NEW sample combined heat maps: "<< outC <<"\n";

        delete h2W32;
        delete h2Best;
    }

    //----------------------------------------------------------------------
    // 4) The EXTRA STEP => 1D red + 2D heat map highlight of top-3 spikes
    //----------------------------------------------------------------------
    {
        // First we build a 1D histogram for “new bestW offset=0”, plus
        // keep track of which IB contributed to each bin => newBestContributors
        TH1F* histNewBest= new TH1F("histNewBest","NEW bestW offset=0;ADC;Count",128,0,16384);
        histNewBest->SetLineColor(kRed);
        histNewBest->SetLineWidth(3);

        // binIndex => vector of (s,b)
        std::map<int, std::vector< std::pair<int,int> >> newBestContributors;

        for (auto &itB : newFitMap) {
            int s= itB.first.first;
            int b= itB.first.second;
            const BoardFitData &bfd = itB.second;

            int bestW= bfd.bestW;
            if(bestW<1) continue; // invalid
            // find bfd.data[ bestW ]
            auto itW= bfd.data.find(bestW);
            if(itW == bfd.data.end()) continue;

            for (auto &xy : itW->second) {
                if(std::fabs(xy.x)<1e-9 && xy.y>0.) {
                    double adcVal= xy.y;
                    int binIdx= histNewBest->FindBin(adcVal);
                    histNewBest->Fill(adcVal);
                    // record contributor
                    newBestContributors[ binIdx ].push_back({s,b});
                }
            }
        }

        // Now identify top-3 spikes in [2000..4000]
        int firstBin= histNewBest->FindBin(2000.);
        int lastBin = histNewBest->FindBin(4000.);
        std::vector< std::pair<int,int> > countsBins; // (count, binIdx)
        for(int ib= firstBin; ib<= lastBin; ib++){
            int c= (int) histNewBest->GetBinContent(ib);
            if(c>0) countsBins.push_back({c, ib});
        }
        // sort descending by c
        std::sort(countsBins.begin(), countsBins.end(),
                  [](auto &a, auto &b){ return a.first> b.first; });

        std::vector< std::tuple<int,int,double> > topBins; // (count,bin,center)
        for(int i=0; i<3 && i<(int)countsBins.size(); i++){
            int c= countsBins[i].first;
            int bn= countsBins[i].second;
            double xc= histNewBest->GetBinCenter(bn);
            topBins.push_back({c,bn,xc});
        }

        // gather all IBs from those top bins
        std::vector<std::pair<int,int>> highlightIBs;
        for (auto &tb : topBins) {
            int binIdx= std::get<1>(tb);
            auto &vec= newBestContributors[ binIdx ];
            highlightIBs.insert(highlightIBs.end(), vec.begin(), vec.end());
        }
        // remove duplicates
        std::sort(highlightIBs.begin(), highlightIBs.end());
        highlightIBs.erase(std::unique(highlightIBs.begin(),highlightIBs.end()),
                           highlightIBs.end());

        // Now the 2×1 canvas => top= 1D red, bottom= 2D color map
        TCanvas cX("cExtra","New BEST: 1D top, 2D bottom", 1200,1200);
        cX.Divide(1,2);

        //-------------------------
        // top => 1D
        //-------------------------
        cX.cd(1);
        gPad->SetRightMargin(0.05); // keep stats box

        double yMax= 1.2* histNewBest->GetMaximum();
        histNewBest->SetMaximum(yMax);
        histNewBest->SetTitle("NEW Sample: best TP (#Delta=0) => 1D distribution");
        histNewBest->GetXaxis()->SetTitle("Average ADC");
        histNewBest->GetYaxis()->SetTitle("Count of IBs");
        histNewBest->Draw("HIST");

        // draw green lines at [2000..4000]
        TLine L1(2000,0,2000,yMax), L2(4000,0,4000,yMax);
        for(TLine* L : {&L1,&L2}){
            L->SetLineColor(kGreen+2);
            L->SetLineStyle(2);
            L->SetLineWidth(2);
            L->Draw("same");
        }

        // “dummy” black line for top-3 peak centers in legend
        TLine peakDummy(0,0,0,0);
        peakDummy.SetLineColor(kBlack);
        peakDummy.SetLineStyle(2);
        peakDummy.SetLineWidth(2);

        // actual black lines for each top bin
        for (auto &tb : topBins) {
            double xC= std::get<2>(tb);
            TLine *vLine= new TLine(xC,0., xC,0.9*yMax);
            vLine->SetLineColor(kBlack);
            vLine->SetLineStyle(2);
            vLine->SetLineWidth(2);
            vLine->Draw("same");
        }

        // legend below the stats box
        TLegend *leg= new TLegend(0.63, 0.62, 0.88, 0.75);
        leg->SetBorderSize(1);
        leg->SetFillColorAlpha(kWhite,0.6);
        leg->SetTextSize(0.03);
        leg->AddEntry(histNewBest, "NEW best TP(#Delta=0)","l");
        leg->AddEntry(&L1, "Target range: 2000..4000","l");
        leg->AddEntry(&peakDummy,"Top-3 peak centers","l");
        leg->Draw();

        // small text block for each top bin
        {
            TLatex tx;
            tx.SetNDC(true);
            tx.SetTextAlign(13); // left-justify
            tx.SetTextSize(0.03);

            double xT=0.63, yT=0.59, dy=0.045;
            for(size_t i=0; i<topBins.size(); i++){
                int count= std::get<0>(topBins[i]);
                double xC= std::get<2>(topBins[i]);
                // e.g. "ADC=2850 => 32 IBs"
                tx.DrawLatex(xT, yT - i*dy,
                             Form("ADC=%.0f => %d IBs", xC, count));
            }
        }

        //-------------------------
        // bottom => 2D
        //-------------------------
        cX.cd(2);
        gPad->SetRightMargin(0.20);

        TH2F* h2= new TH2F("h2","NEW bestW offset=0 (2D);#phi index;#eta index",
                           256,0,256, 96,0,96);
        h2->SetStats(false);
        h2->GetZaxis()->SetTitle("avg ADC");

        // fill from getBestOffset0(false => new)
        for (auto &itB : newFitMap) {
            int s= itB.first.first;
            int b= itB.first.second;
            double val= getBestOffset0(false,s,b);
            if(val<0) continue;

            int localSec= (s % 32);
            int phiMin= 8*localSec, phiMax= phiMin+7;
            int etaMin=0, etaMax=0;
            if(s<32) {
                etaMin=48 + 8*b;
                etaMax= etaMin+7;
            } else {
                etaMin=40 - 8*b;
                etaMax= etaMin+7;
            }
            if(etaMin>etaMax) std::swap(etaMin,etaMax);

            for(int xx= phiMin; xx<=phiMax; xx++){
                for(int yy= etaMin; yy<=etaMax; yy++){
                    h2->SetBinContent(xx+1,yy+1, val);
                }
            }
        }
        h2->Draw("COLZ");

        // label each IB => "s.b"
        for (auto &itB : newFitMap) {
            int s= itB.first.first;
            int b= itB.first.second;
            double val= getBestOffset0(false,s,b);
            if(val<0) continue;

            int localSec= (s % 32);
            int phiMin= 8*localSec, phiMax= phiMin+7;
            int etaMin=0, etaMax=0;
            if(s<32){
                etaMin=48+8*b; etaMax=etaMin+7;
            } else {
                etaMin=40- 8*b; etaMax= etaMin+7;
            }
            if(etaMin>etaMax) std::swap(etaMin,etaMax);

            double cx= 0.5*(phiMin + phiMax +1);
            double cy= 0.5*(etaMin + etaMax +1);

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
                etaMin=48+8*b; etaMax= etaMin+8;
            } else {
                etaMin=40-8*b; etaMax= etaMin+8;
            }
            if(etaMin>etaMax) std::swap(etaMin,etaMax);

            double cx= 0.5*(phiMin + phiMax);
            double cy= 0.5*(etaMin + etaMax);

            TMarker *mk= new TMarker(cx,cy,24); // open circle
            mk->SetMarkerColor(kRed);
            mk->SetMarkerSize(2.5);
            mk->Draw();
        }

        // final
        std::string outEx= summary + "/newBest_1D_top_2D_bottom_highlight.png";
        cX.SaveAs(outEx.c_str());
        std::cout<<"[INFO] => Created extra figure => "<< outEx <<"\n";

        delete h2;
        delete histNewBest;
    }

    //----------------------------------------------------------------------
    // 5) Print the tallies => old_total, new_total, etc.
    //----------------------------------------------------------------------
    int old_total= old_below + old_inrange + old_above + old_none;
    int new_total= new_below + new_inrange + new_above + new_none;

    std::cout << "\n[DEBUG] OLD(32,0) color-coded map => total IB= "<< old_total << "\n"
              << "  below=  " << old_below <<" => "<< (100.*old_below/old_total)  <<"%\n"
              << "  inRange=" << old_inrange<<" => "<< (100.*old_inrange/old_total)<<"%\n"
              << "  above=  " << old_above  <<" => "<< (100.*old_above/old_total)  <<"%\n"
              << "  none=   " << old_none   <<" => "<< (100.*old_none/old_total)   <<"%\n";

    std::cout << "[DEBUG] NEW(32,0) color-coded map => total IB= "<< new_total << "\n"
              << "  below=  " << new_below <<" => "<< (100.*new_below/new_total)  <<"%\n"
              << "  inRange=" << new_inrange<<" => "<< (100.*new_inrange/new_total)<<"%\n"
              << "  above=  " << new_above  <<" => "<< (100.*new_above/new_total)  <<"%\n"
              << "  none=   " << new_none   <<" => "<< (100.*new_none/new_total)   <<"%\n";
}


void outOfSigmaAmplitudeCheck(
    const std::vector<Kval> &kRatios,     // has kOld,kNew for each board
    double meanN,                         // the global mean(kNew) from your fits
    double sigmaN,                        // the global sigma(kNew) from your fits
    const std::map< std::pair<int,int>, BoardFitData > &newFitMap,
    const std::map< std::pair<int,int>, BoardFitData > &oldFitMap,
    const std::string &outputPNG,         // e.g. "/path/amplitudeCheck_outliers_twopanels.png"
    const std::string &extraLogDir        // for printing usage summary in logs
)
{
    //-----------------------------------------------------------------------
    // 1) Classify boards into in-sigma vs. out-of-sigma => gather amplitudes
    //-----------------------------------------------------------------------
    std::vector<double> amps_inSigma;
    std::vector<double> amps_outSigma;

    std::cout << "\n[CHECK] => Checking amplitude variation for 'out-of-sigma' boards...\n";
    for (auto &kk : kRatios)
    {
        // skip if kNew <= 0 => invalid
        if (kk.kNew <= 0.) continue;

        // check how far kNew is from meanN
        double diff = std::fabs(kk.kNew - meanN);
        bool isOutlier = (diff > sigmaN);

        // Read amplitude from newFitMap
        auto itB = newFitMap.find({kk.sector, kk.ib});
        if (itB == newFitMap.end()) {
            // not found => skip
            continue;
        }

        double aN = itB->second.amplitude; // "amplitude" from bestW fit
        if (aN <= 0.) {
            // amplitude <=0 => skip
            continue;
        }

        // classify amplitude
        if (isOutlier) {
            amps_outSigma.push_back(aN);
            std::cout << "  [OUTLIER] sector=" << kk.sector
                      << ", IB=" << kk.ib
                      << ", kNew=" << kk.kNew
                      << ", amplitude=" << aN
                      << std::endl;
        }
        else {
            amps_inSigma.push_back(aN);
        }
    }

    //-----------------------------------------------------------------------
    // 2) Build histograms => hA_in, hA_out
    //-----------------------------------------------------------------------
    double maxA = 0.;
    for (auto v : amps_inSigma)  { if (v>maxA) maxA=v; }
    for (auto v : amps_outSigma) { if (v>maxA) maxA=v; }
    if (maxA < 1e-7) maxA = 1e-5;

    TH1F hA_in ("hA_in","Intercept for in-sigma kNew;Intercept (amplitude);Count",
                50, 0., 1.2*maxA);
    TH1F hA_out("hA_out","Intercept for out-of-sigma kNew;Intercept (amplitude);Count",
                50, 0., 1.2*maxA);

    for (auto v : amps_inSigma)  hA_in.Fill(v);
    for (auto v : amps_outSigma) hA_out.Fill(v);

    std::cout << "\n[CHECK] => # in-sigma (kNew) boards  = "
              << amps_inSigma.size()  << std::endl
              << "          # out-of-sigma (kNew) boards = "
              << amps_outSigma.size() << std::endl;

    //-----------------------------------------------------------------------
    // 3) Compute arithmetic mean/median for each group
    //-----------------------------------------------------------------------
    double meanIn=0., medianIn=0.;
    if (!amps_inSigma.empty()) {
        meanIn   = TMath::Mean(   (int)amps_inSigma.size(),  &amps_inSigma[0] );
        medianIn = TMath::Median( (int)amps_inSigma.size(), &amps_inSigma[0] );
    }

    double meanOut=0., medianOut=0.;
    if (!amps_outSigma.empty()) {
        meanOut   = TMath::Mean(   (int)amps_outSigma.size(),  &amps_outSigma[0] );
        medianOut = TMath::Median( (int)amps_outSigma.size(), &amps_outSigma[0] );
    }

    //-----------------------------------------------------------------------
    // 4) Two-panel canvas => left=in-sigma, right=out-of-sigma
    //-----------------------------------------------------------------------
    TCanvas cCheck("cCheck","Amplitude Variation for Out-of-Sigma kNew",1200,600);
    cCheck.Divide(2,1);

    // left => in-sigma
    cCheck.cd(1);
    hA_in.SetLineColor(kBlue);
    hA_in.SetLineWidth(2);
    hA_in.Draw("HIST");
    {
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextSize(0.03);
        double xPos= 0.55, yPos= 0.78;

        lat.DrawLatex(xPos, yPos, Form("k_{#mu} = %.3g", meanN));   yPos-= 0.04;
        lat.DrawLatex(xPos, yPos, Form("k_{#sigma} = %.3g", sigmaN));yPos-= 0.055;

        lat.DrawLatex(xPos, yPos, Form("#IBs in-sigma = %d", (int)amps_inSigma.size()));
        yPos-= 0.04;
        lat.DrawLatex(xPos, yPos, Form("Arithmetic Mean=%.1f", meanIn));
        yPos-= 0.04;
        lat.DrawLatex(xPos, yPos, Form("Median=%.1f", medianIn));
    }

    // right => out-of-sigma
    cCheck.cd(2);
    hA_out.SetLineColor(kRed);
    hA_out.SetLineWidth(2);
    hA_out.Draw("HIST");
    {
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextSize(0.03);
        double xPos= 0.55, yPos= 0.78;

        lat.DrawLatex(xPos, yPos, Form("k_{#mu} = %.3g", meanN));   yPos-= 0.04;
        lat.DrawLatex(xPos, yPos, Form("k_{#sigma} = %.3g", sigmaN));yPos-= 0.055;

        lat.DrawLatex(xPos, yPos, Form("#IBs out-of-sigma = %d", (int)amps_outSigma.size()));
        yPos-= 0.04;
        lat.DrawLatex(xPos, yPos, Form("Arithmetic Mean=%.1f", meanOut));
        yPos-= 0.04;
        lat.DrawLatex(xPos, yPos, Form("Median=%.1f", medianOut));
    }

    // Save the canvas
    cCheck.SaveAs(outputPNG.c_str());
    std::cout << "[CHECK] => Created plot: " << outputPNG << std::endl;
    std::cout << "[CHECK] => Done with out-of-sigma amplitude cross-check!\n\n";

    //-----------------------------------------------------------------------
    // 5) Best pulse width usage summary in [24..40] for old/new bestW
    //-----------------------------------------------------------------------
    {
        // 1) Prepare counters
        std::map<int,int> countOld;
        std::map<int,int> countNew;
        for(int w=24; w<=40; w++){
            countOld[w] = 0;
            countNew[w] = 0;
        }

        // 2) Loop oldFitMap => old bestW
        for (auto &pr : oldFitMap) {
            int wUsed = pr.second.bestW;
            if (wUsed>=24 && wUsed<=40) {
                countOld[wUsed] += 1;
            }
        }

        // 3) Loop newFitMap => new bestW
        for (auto &pr : newFitMap) {
            int wUsed = pr.second.bestW;
            if (wUsed>=24 && wUsed<=40) {
                countNew[wUsed] += 1;
            }
        }

        // 4) ANSI color codes
        const char* RESET  = "\033[0m";
        const char* GREEN  = "\033[32m";
        const char* RED    = "\033[31m";
        const char* CYAN   = "\033[36m";
        const char* YELLOW = "\033[93m";

        // Print heading
        std::cout << "\n" << YELLOW
                  << "=== BEST PULSE WIDTH USAGE SUMMARY (24..40 ns) ==="
                  << RESET << "\n"
                  << "   (How many IBs used each width as 'best'?)\n\n";

        std::cout << CYAN
                  << "  Width(ns)   |  # OLD   |  # NEW  "
                  << RESET << "\n"
                  << "  -----------   --------   -------\n";

        // 5) Print rows
        for(int w=24; w<=40; w++){
            int nOld= countOld[w];
            int nNew= countNew[w];

            std::string colO= (nOld>0 ? GREEN : RED);
            std::string colN= (nNew>0 ? GREEN : RED);

            char line[256];
            std::snprintf(line,sizeof(line),
                          "     %2d       |   %s%4d%s   |   %s%4d%s",
                          w,
                          colO.c_str(), nOld, RESET,
                          colN.c_str(), nNew, RESET);
            std::cout << line << "\n";
        }

        std::cout << "\n" << YELLOW
                  << "=== END of BEST PULSE WIDTH USAGE SUMMARY ==="
                  << RESET << "\n\n";
    }
}


void doAnalysis(
    const std::string &oldRunCSV,
    const std::string &oldInputDir,
    const std::string &newRunCSV,
    const std::string &newInputDir,
    const std::string &vopCSV,
    const std::string &outDir
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
              << " vopCSV="      << vopCSV       << "\n"
              << " outDir="      << outDir       << "\n\n";

    if (gSystem->AccessPathName(outDir.c_str())) {
        gSystem->Exec(("mkdir -p " + outDir).c_str());
    }
    std::string multiW  = outDir + "/multiW";
    std::string summary = outDir + "/summary";
    if (gSystem->AccessPathName(multiW.c_str())) {
        gSystem->Exec(("mkdir -p " + multiW).c_str());
    }
    if (gSystem->AccessPathName(summary.c_str())) {
        gSystem->Exec(("mkdir -p " + summary).c_str());
    }

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

    produceCompare32vsBest(dataMap, bestWmap, summary);

    produceKAndAmplitudePlots(oldFitMap, newFitMap, summary);
    
    produceColorCodedMapsAndExtra(
        oldFitMap,
        newFitMap,
        sectorIBmap,
        summary  // outDir + "/summary"
    );

    std::string outPng = "/Users/patsfan753/Desktop/PMTgainsAna/output/amplitudeCheck_outliers_twopanels.png";

    // Then just call:
    outOfSigmaAmplitudeCheck(
        kRatios,
        meanN,
        sigmaN,
        newFitMap,
        oldFitMap,
        outPng,
        summary  // or another path if you prefer
    );
    
    produceNewPerIBOverlays(
        newFitMap,        // pass map of (sector,ib)->BoardFitData
        amplitudeRatios,  // pass existing vector of amplitude info
        sectorIBmap,      // so we know which (sector,IB) to iterate
        "/Users/patsfan753/Desktop/PMTgainsAna/output/NewDataGainCurves"  // output folder
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
        basePath + "vop.csv",
        "/Users/patsfan753/Desktop/PMTgainsAna/output/compareOldNew"
    );
 
  // End with the original final message:
  std::cout<<"[INFO] ---------- ANALYSIS COMPLETE ---------.\n";
}
