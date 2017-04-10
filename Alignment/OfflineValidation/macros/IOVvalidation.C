#include <vector>
#include <string>
#include <iostream>
#include "TCanvas.h"
#include "TH1D.h"
#include <TSystem.h>
#include <sstream>
#include "TLegend.h"
#include <TStyle.h>
#include <map>
#include <algorithm>
#include <TLine.h>
#include <TFrame.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/foreach.hpp>
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TBranch.h"
#include "TLeaf.h"

namespace pt = boost::property_tree;

class timeplot
{
public:
    struct IOVData
    {
        std::string IOVlabel;
        std::string filename;
        std::string treename;
        double rms[8], mu[8], deltamu[8], rmserror[8], muerror[8], deltamuerror[8];
        int linecolor;
    };
    struct Changeline
    {
        std::string changename;
        std::string previousrun;
        int linecolor2;
    };
    std::string subdetectors="";
    std::string drawparameter="";
    void draw(const std::string& subdet, const std::string& parameter = "rms,mu,deltamu");
    void setoutputdir(const std::string& dir = ".");
    void addIOV(const std::string& IOVid, const std::string& filepath, const std::string& tree, const std::string& alignment, const int colorline);
    void sequence();
    void addlineofchange(const std::string& nameofchange, const std::string& previousrun, const int colorline);
private:
    std::string outputdir;
    std::map<std::string, std::vector<IOVData>> alignments;
    std::map<std::string, int> seq;
    std::vector<Changeline> lines;
    std::map<std::string, int> sdet = {{"BPIX",0},{"FPIX",1},{"TIB",2},{"TID",3},{"TOB",4},{"TEC",5},{"BPIX_Y",6},{"FPIX_Y",7}};
};

void timeplot::draw(const std::string& subdet, const std::string& parameter)
{
    size_t multipara=parameter.find(",");
    if(multipara!=std::string::npos)
    {
        std::string substring1 = parameter.substr(0,multipara);
        std::string substring2 = parameter.substr(multipara+1,std::string::npos);
        draw(subdet, substring1);
        draw(subdet, substring2);
        return;
    }
    size_t multidet=subdet.find(",");
    if(multidet!=std::string::npos)
    {
        std::string substring3 = subdet.substr(0,multidet);
        std::string substring4 = subdet.substr(multidet+1,std::string::npos);
        draw(substring3, parameter);
        draw(substring4, parameter);
        return;
    }
    TCanvas c("canv", "canv");
    gStyle->SetOptStat(0);
    int numberIOV=static_cast<int>(seq.size());
    std::vector<TH1D*> rmsvec;
    std::vector<TH1D*> muvec;
    std::vector<TH1D*> deltamuvec;
    std::ostringstream plotName;
    plotName << outputdir << "/";
    TLegend* leg=new TLegend(0.6,0.85,0.8,1.0,"Alignments");
    TLegend* leg2=new TLegend(0.8,0.85,1.0,1.0,"Detector changes");
    leg->SetTextSize(0.025);
    leg2->SetTextSize(0.025);
    
    std::vector<int> sorted;
    int a=1;
    for(std::map<std::string, int>::iterator itseq=seq.begin(); itseq!=seq.end(); itseq++)
    {
        sorted.push_back(std::stoi(itseq->first));
        a++;
    }
    std::sort(sorted.begin(),sorted.end());
    double ymax[3]={0,0,0};
    double ymin[3]={0,0,0};
    for(std::map<std::string, std::vector<IOVData>>::iterator q=alignments.begin(); q!=alignments.end(); q++)
    {
        std::vector<IOVData> w=alignments[q->first];
        for(std::vector<IOVData>::iterator q2=alignments[q->first].begin(); q2!=alignments[q->first].end(); q2++)
        {
            int e=q2-alignments[q->first].begin();
            if((w[e].rms[sdet[subdet]]+w[e].rmserror[sdet[subdet]])>ymax[0])
            {
                ymax[0]=w[e].rms[sdet[subdet]]+w[e].rmserror[sdet[subdet]];
            }
            if((w[e].mu[sdet[subdet]]+w[e].muerror[sdet[subdet]])>ymax[1])
            {
                ymax[1]=w[e].mu[sdet[subdet]]+w[e].muerror[sdet[subdet]];
            }
            if((w[e].deltamu[sdet[subdet]]+w[e].deltamuerror[sdet[subdet]])>ymax[2])
            {
                ymax[2]=w[e].deltamu[sdet[subdet]]+w[e].deltamuerror[sdet[subdet]];
            }
            if((w[e].rms[sdet[subdet]]-w[e].rmserror[sdet[subdet]])<ymin[0])
            {
                ymin[0]=w[e].rms[sdet[subdet]]-w[e].rmserror[sdet[subdet]];
            }
            if((w[e].mu[sdet[subdet]]-w[e].muerror[sdet[subdet]])<ymin[1])
            {
                ymin[1]=w[e].mu[sdet[subdet]]-w[e].muerror[sdet[subdet]];
            }
            if((w[e].deltamu[sdet[subdet]]-w[e].deltamuerror[sdet[subdet]])<ymin[2])
            {
                ymin[2]=w[e].deltamu[sdet[subdet]]-w[e].deltamuerror[sdet[subdet]];
            }
        }
    }
    for(std::map<std::string, std::vector<IOVData>>::iterator it=alignments.begin(); it!=alignments.end(); it++)
    {
        int bin;
        std::vector<IOVData> ali=alignments[it->first];
        if(parameter=="rms")
        {
	    TUUID id1;
            TH1D* rms=new TH1D(id1.AsString(), subdet.c_str(), numberIOV, 0, numberIOV);
            rms->SetLineColor(ali[0].linecolor);
            for(std::vector<IOVData>::iterator it2=alignments[it->first].begin(); it2!=alignments[it->first].end(); it2++)
            {
                int i=it2-alignments[it->first].begin();
                bin=seq[ali[i].IOVlabel]+1;
                rms->SetBinContent(bin, ali[i].rms[sdet[subdet]]);
                rms->SetBinError(bin, ali[i].rmserror[sdet[subdet]]);
            }
            for(int b=0; b<numberIOV; b++)
            {
                rms->GetXaxis()->SetBinLabel(b+1, to_string(sorted[b]).c_str());
            }
            rmsvec.push_back(rms);
            leg->AddEntry(rms, it->first.c_str(), "L");
        }
        if(parameter=="mu")
        {
	    TUUID id2;
            TH1D* mu=new TH1D(id2.AsString(), subdet.c_str(), numberIOV, 0, numberIOV);
            mu->SetLineColor(ali[0].linecolor);
            for(std::vector<IOVData>::iterator it2=alignments[it->first].begin(); it2!=alignments[it->first].end(); it2++)
            {
                int i=it2-alignments[it->first].begin();
                bin=seq[ali[i].IOVlabel]+1;
                mu->SetBinContent(bin, ali[i].mu[sdet[subdet]]);
                mu->SetBinError(bin, ali[i].muerror[sdet[subdet]]);
            }
            for(int b=0; b<numberIOV; b++)
            {
                mu->GetXaxis()->SetBinLabel(b+1, to_string(sorted[b]).c_str());
            }
            muvec.push_back(mu);
            leg->AddEntry(mu, it->first.c_str(), "L");
        }
        if(parameter=="deltamu")
        {
	    TUUID id3;
            TH1D* deltamu=new TH1D(id3.AsString(), subdet.c_str(), numberIOV, 0, numberIOV);
            deltamu->SetLineColor(ali[0].linecolor);
            for(std::vector<IOVData>::iterator it2=alignments[it->first].begin(); it2!=alignments[it->first].end(); it2++)
            {
                int i=it2-alignments[it->first].begin();
                bin=seq[ali[i].IOVlabel]+1;
                deltamu->SetBinContent(bin, ali[i].deltamu[sdet[subdet]]);
                deltamu->SetBinError(bin, ali[i].deltamuerror[sdet[subdet]]);
            }
            for(int b=0; b<numberIOV; b++)
            {
                deltamu->GetXaxis()->SetBinLabel(b+1, to_string(sorted[b]).c_str());
            }
            deltamuvec.push_back(deltamu);
            leg->AddEntry(deltamu, it->first.c_str(), "L");
        }
    }
    std::string linetyps="";
    int linestyle=1;
    if(parameter=="rms")
    {
        rmsvec[0]->SetMaximum(ymax[0]*1.5);
        rmsvec[0]->SetMinimum(ymin[0]*1.2);
        rmsvec[0]->Draw("E");
        c.SetBottomMargin(0.15);
        rmsvec[0]->GetXaxis()->SetTitle("IOV number");
        rmsvec[0]->GetYaxis()->SetTitle("rms in #mum");
        rmsvec[0]->GetXaxis()->SetTitleOffset(1.8);
	rmsvec[0]->LabelsOption("v");
        for(std::vector<TH1D*>::iterator it3=rmsvec.begin()+1; it3!=rmsvec.end(); it3++)
        {
            int i2=it3-rmsvec.begin();
            rmsvec[i2]->Draw("SAME E");
	    c.SetBottomMargin(0.15);
        }
        c.Update();
        for(std::vector<Changeline>::iterator lineit=lines.begin(); lineit!=lines.end(); lineit++)
        {
            int lit=lineit-lines.begin();
            int binnumber=seq[lines[lit].previousrun];
            int xline=(rmsvec[0]->GetBinCenter(binnumber+2)+rmsvec[0]->GetBinCenter(binnumber+1))/2.;
            TLine *line=new TLine(xline,c.GetFrame()->GetY1(),xline,c.GetFrame()->GetY2());
            line->SetLineColor(lines[lit].linecolor2);
	    if(linetyps.find(lines[lit].changename.c_str())==std::string::npos)
	    {
	      leg2->AddEntry(line, lines[lit].changename.c_str(), "L");
	      linetyps.append(lines[lit].changename);
	      linestyle+=1;
	    }
	    line->SetLineStyle(linestyle);
            line->Draw();
        }
        leg->Draw();
        leg2->Draw();
        gPad->RedrawAxis();
        c.SetBottomMargin(0.15);
        c.Print((plotName.str()+subdet+"_rms.pdf").c_str());
    }
    if(parameter=="mu")
    {
        muvec[0]->SetMaximum(ymax[1]*1.5);
        muvec[0]->SetMinimum(ymin[1]*1.2);
        muvec[0]->Draw("E");
        c.SetBottomMargin(0.15);
        muvec[0]->GetXaxis()->SetTitle("IOV number");
        muvec[0]->GetYaxis()->SetTitle("#mu in #mum");
        muvec[0]->GetXaxis()->SetTitleOffset(1.8);
	muvec[0]->LabelsOption("v");
        for(std::vector<TH1D*>::iterator it4=muvec.begin()+1; it4!=muvec.end(); it4++)
        {
            int i3=it4-muvec.begin();
            muvec[i3]->Draw("SAME E");
	    c.SetBottomMargin(0.15);
        }
        c.Update();
        for(std::vector<Changeline>::iterator lineit=lines.begin(); lineit!=lines.end(); lineit++)
        {
            int lit=lineit-lines.begin();
            int binnumber=seq[lines[lit].previousrun];
            int xline=(muvec[0]->GetBinCenter(binnumber+2)+muvec[0]->GetBinCenter(binnumber+1))/2.;
            TLine *line=new TLine(xline,c.GetFrame()->GetY1(),xline,c.GetFrame()->GetY2());
            line->SetLineColor(lines[lit].linecolor2);
	    if(linetyps.find(lines[lit].changename.c_str())==std::string::npos)
	    {
	      leg2->AddEntry(line, lines[lit].changename.c_str(), "L");
	      linetyps.append(lines[lit].changename);
	      linestyle+=1;
	    }
	    line->SetLineStyle(linestyle);
            line->Draw();
        }
        leg->Draw();
        leg2->Draw();
        gPad->RedrawAxis();
        c.SetBottomMargin(0.15);
        c.Print((plotName.str()+subdet+"_mu.pdf").c_str());
    }
    if(parameter=="deltamu")
    {
        deltamuvec[0]->SetMaximum(ymax[2]*1.5);
        deltamuvec[0]->SetMinimum(ymin[2]*1.2);
        deltamuvec[0]->Draw("E");
        c.SetBottomMargin(0.15);
        deltamuvec[0]->GetXaxis()->SetTitle("IOV number");
        deltamuvec[0]->GetYaxis()->SetTitle("#Delta#mu in #mum");
        deltamuvec[0]->GetXaxis()->SetTitleOffset(1.8);
	deltamuvec[0]->LabelsOption("v");
        for(std::vector<TH1D*>::iterator it5=deltamuvec.begin()+1; it5!=deltamuvec.end(); it5++)
        {
            int i4=it5-deltamuvec.begin();
            deltamuvec[i4]->Draw("SAME E");
	    c.SetBottomMargin(0.15);
        }
        c.Update();
        for(std::vector<Changeline>::iterator lineit=lines.begin(); lineit!=lines.end(); lineit++)
        {
            int lit=lineit-lines.begin();
            int binnumber=seq[lines[lit].previousrun];
            int xline=(deltamuvec[0]->GetBinCenter(binnumber+2)+deltamuvec[0]->GetBinCenter(binnumber+1))/2.;
            TLine *line=new TLine(xline,c.GetFrame()->GetY1(),xline,c.GetFrame()->GetY2());
            line->SetLineColor(lines[lit].linecolor2);
	    line->SetLineStyle(2);
            line->Draw();
	    if(linetyps.find(lines[lit].changename.c_str())==std::string::npos)
	    {
	      leg2->AddEntry(line, lines[lit].changename.c_str(), "L");
	      linetyps.append(lines[lit].changename);
	      linestyle+=1;
	    }
	    line->SetLineStyle(linestyle);
            line->Draw();
        }
        leg->Draw();
        leg2->Draw();
        gPad->RedrawAxis();
        c.SetBottomMargin(0.15);
        c.Print((plotName.str()+subdet+"_deltamu.pdf").c_str());
    }
}

void timeplot::setoutputdir(const std::string& dir )
{
    std::cout <<"'"<< outputdir <<"' = "<< dir << std::endl;
    outputdir = dir;
    gSystem->mkdir(outputdir.data(), true);
}

void timeplot::addIOV(const std::string& IOVid, const std::string&  filepath,  const std::string& tree, const std::string& alignment, const int colorline)
{
    IOVData newIOV;
    newIOV.IOVlabel=IOVid;
    newIOV.filename=filepath;
    newIOV.treename=tree;
    TFile *para=new TFile(newIOV.filename.c_str(),"read");
    TTree *t1=(TTree*)para->Get(newIOV.treename.c_str());
    t1->GetEntry(0);
    newIOV.rms[0]=t1->GetBranch("BPIX")->GetLeaf("rms")->GetValue(0);
    newIOV.rms[1]=t1->GetBranch("FPIX")->GetLeaf("rms")->GetValue(0);
    newIOV.rms[2]=t1->GetBranch("TIB")->GetLeaf("rms")->GetValue(0);
    newIOV.rms[3]=t1->GetBranch("TID")->GetLeaf("rms")->GetValue(0);
    newIOV.rms[4]=t1->GetBranch("TOB")->GetLeaf("rms")->GetValue(0);
    newIOV.rms[5]=t1->GetBranch("TEC")->GetLeaf("rms")->GetValue(0);
    newIOV.rms[6]=t1->GetBranch("BPIX_Y")->GetLeaf("rms")->GetValue(0);
    newIOV.rms[7]=t1->GetBranch("FPIX_Y")->GetLeaf("rms")->GetValue(0);
    newIOV.mu[0]=t1->GetBranch("BPIX")->GetLeaf("mu")->GetValue(0);
    newIOV.mu[1]=t1->GetBranch("FPIX")->GetLeaf("mu")->GetValue(0);
    newIOV.mu[2]=t1->GetBranch("TIB")->GetLeaf("mu")->GetValue(0);
    newIOV.mu[3]=t1->GetBranch("TID")->GetLeaf("mu")->GetValue(0);
    newIOV.mu[4]=t1->GetBranch("TOB")->GetLeaf("mu")->GetValue(0);
    newIOV.mu[5]=t1->GetBranch("TEC")->GetLeaf("mu")->GetValue(0);
    newIOV.mu[6]=t1->GetBranch("BPIX_Y")->GetLeaf("mu")->GetValue(0);
    newIOV.mu[7]=t1->GetBranch("FPIX_Y")->GetLeaf("mu")->GetValue(0);
    newIOV.deltamu[0]=t1->GetBranch("BPIX")->GetLeaf("deltamu")->GetValue(0);
    newIOV.deltamu[1]=t1->GetBranch("FPIX")->GetLeaf("deltamu")->GetValue(0);
    newIOV.deltamu[2]=t1->GetBranch("TIB")->GetLeaf("deltamu")->GetValue(0);
    newIOV.deltamu[3]=0;
    newIOV.deltamu[4]=t1->GetBranch("TOB")->GetLeaf("deltamu")->GetValue(0);
    newIOV.deltamu[5]=0;
    newIOV.deltamu[6]=t1->GetBranch("BPIX_Y")->GetLeaf("deltamu")->GetValue(0);
    newIOV.deltamu[7]=t1->GetBranch("FPIX_Y")->GetLeaf("deltamu")->GetValue(0);
    newIOV.rmserror[0]=t1->GetBranch("BPIX")->GetLeaf("rmserror")->GetValue(0);
    newIOV.rmserror[1]=t1->GetBranch("FPIX")->GetLeaf("rmserror")->GetValue(0);
    newIOV.rmserror[2]=t1->GetBranch("TIB")->GetLeaf("rmserror")->GetValue(0);
    newIOV.rmserror[3]=t1->GetBranch("TID")->GetLeaf("rmserror")->GetValue(0);
    newIOV.rmserror[4]=t1->GetBranch("TOB")->GetLeaf("rmserror")->GetValue(0);
    newIOV.rmserror[5]=t1->GetBranch("TEC")->GetLeaf("rmserror")->GetValue(0);
    newIOV.rmserror[6]=t1->GetBranch("BPIX_Y")->GetLeaf("rmserror")->GetValue(0);
    newIOV.rmserror[7]=t1->GetBranch("FPIX_Y")->GetLeaf("rmserror")->GetValue(0);
    newIOV.muerror[0]=t1->GetBranch("BPIX")->GetLeaf("muerror")->GetValue(0);
    newIOV.muerror[1]=t1->GetBranch("FPIX")->GetLeaf("muerror")->GetValue(0);
    newIOV.muerror[2]=t1->GetBranch("TIB")->GetLeaf("muerror")->GetValue(0);
    newIOV.muerror[3]=t1->GetBranch("TID")->GetLeaf("muerror")->GetValue(0);
    newIOV.muerror[4]=t1->GetBranch("TOB")->GetLeaf("muerror")->GetValue(0);
    newIOV.muerror[5]=t1->GetBranch("TEC")->GetLeaf("muerror")->GetValue(0);
    newIOV.muerror[6]=t1->GetBranch("BPIX_Y")->GetLeaf("muerror")->GetValue(0);
    newIOV.muerror[7]=t1->GetBranch("FPIX_Y")->GetLeaf("muerror")->GetValue(0);
    newIOV.deltamuerror[0]=t1->GetBranch("BPIX")->GetLeaf("deltamuerror")->GetValue(0);
    newIOV.deltamuerror[1]=t1->GetBranch("FPIX")->GetLeaf("deltamuerror")->GetValue(0);
    newIOV.deltamuerror[2]=t1->GetBranch("TIB")->GetLeaf("deltamuerror")->GetValue(0);
    newIOV.deltamuerror[3]=0;
    newIOV.deltamuerror[4]=t1->GetBranch("TOB")->GetLeaf("deltamuerror")->GetValue(0);
    newIOV.deltamuerror[5]=0;
    newIOV.deltamuerror[6]=t1->GetBranch("BPIX_Y")->GetLeaf("deltamuerror")->GetValue(0);
    newIOV.deltamuerror[7]=t1->GetBranch("FPIX_Y")->GetLeaf("deltamuerror")->GetValue(0);
    delete para;
    newIOV.linecolor=colorline;
    alignments[alignment].push_back(newIOV);
}

void timeplot::sequence()
{
    std::vector<int> times;
    for(std::map<std::string, std::vector<IOVData>>::iterator it=alignments.begin(); it!=alignments.end(); it++)
    {
        for(std::vector<IOVData>::iterator it2=alignments[it->first].begin(); it2!=alignments[it->first].end(); it2++)
        {
            int i=it2-alignments[it->first].begin();
            std::vector<IOVData> al=alignments[it->first];
            int time=std::stoi(al[i].IOVlabel);
            bool newtime=true;
            for(std::vector<int>::iterator it3=times.begin(); it3!=times.end(); it3++)
            {
                int i2=it3-times.begin();
                if(time==times[i2])
                {
                    newtime=false;
                }
            }
            if(newtime)
            {
                times.push_back(time);
            }
        }
    }
    std::sort(times.begin(), times.end());
    for(std::map<std::string, std::vector<IOVData>>::iterator it4=alignments.begin(); it4!=alignments.end(); it4++)
    {
        for(std::vector<IOVData>::iterator it5=alignments[it4->first].begin(); it5!=alignments[it4->first].end(); it5++)
        {
            int i2=it5-alignments[it4->first].begin();
            std::vector<IOVData> al2=alignments[it4->first];
            for(std::vector<int>::iterator it6=times.begin(); it6!=times.end(); it6++)
            {
                int i3=it6-times.begin();
                if(times[i3]==std::stoi(al2[i2].IOVlabel))
                {
                    seq[al2[i2].IOVlabel]=i3;
                }
            }
        }
    }       
}

void timeplot::addlineofchange(const std::string& nameofchange, const std::string& previousrun, const int colorline)
{
    Changeline newline;
    newline.changename=nameofchange;
    newline.previousrun=previousrun;
    newline.linecolor2=colorline;
    lines.push_back(newline);
}

void IOVvalidation(void)
{
    gROOT->SetBatch(kTRUE);
    pt::ptree tree;
    pt::ini_parser::read_ini("config.ini", tree);
    timeplot p;
    std::string alignmentsToAdd, iovsToAdd, rootpathsToAdd, treenamesToAdd, colorsToAdd, descriptionsToAdd, prevToAdd, colorsToAdd2;
    alignmentsToAdd=tree.get<std::string>("alignment.name");
    iovsToAdd=tree.get<std::string>("alignment.IOV");
    rootpathsToAdd=tree.get<std::string>("alignment.rootpath");
    treenamesToAdd=tree.get<std::string>("alignment.treename");
    colorsToAdd=tree.get<std::string>("alignment.color");
    descriptionsToAdd=tree.get<std::string>("detectorchange.description");
    prevToAdd=tree.get<std::string>("detectorchange.previousrun");
    colorsToAdd2=tree.get<std::string>("detectorchange.color");
    while(alignmentsToAdd.find(",")!=std::string::npos)
    {
        std::size_t pos=alignmentsToAdd.find(",");
        std::string addalingment=alignmentsToAdd.substr(0,pos);
        alignmentsToAdd=alignmentsToAdd.substr(pos+1);
        pos=iovsToAdd.find(",");
        std::string addiovs=iovsToAdd.substr(0,pos);
        iovsToAdd=iovsToAdd.substr(pos+1);
        pos=rootpathsToAdd.find(",");
        std::string addroot=rootpathsToAdd.substr(0,pos);
        rootpathsToAdd=rootpathsToAdd.substr(pos+1);
        pos=treenamesToAdd.find(",");
        std::string addtree=treenamesToAdd.substr(0,pos);
        treenamesToAdd=treenamesToAdd.substr(pos+1);
        pos=colorsToAdd.find(",");
        std::string stringcol=colorsToAdd.substr(0,pos);
        colorsToAdd=colorsToAdd.substr(pos+1);
        int addcolor=std::stoi(stringcol);
        p.addIOV(addiovs, addroot, addtree, addalingment, addcolor);
    }
    int addcolor2=std::stoi(colorsToAdd);
    p.addIOV(iovsToAdd, rootpathsToAdd, treenamesToAdd, alignmentsToAdd, addcolor2);
    while(descriptionsToAdd.find(",")!=std::string::npos)
    {
        std::size_t pos=descriptionsToAdd.find(",");
        std::string adddesc=descriptionsToAdd.substr(0,pos);
        descriptionsToAdd=descriptionsToAdd.substr(pos+1);
        pos=prevToAdd.find(",");
        std::string addprev=prevToAdd.substr(0,pos);
        prevToAdd=prevToAdd.substr(pos+1);
        pos=colorsToAdd2.find(",");
        std::string stringcol2=colorsToAdd2.substr(0,pos);
        colorsToAdd2=colorsToAdd2.substr(pos+1);
        int addcolor=std::stoi(stringcol2);
        p.addlineofchange(adddesc, addprev, addcolor);
    }
    int addcolor3=std::stoi(colorsToAdd2);
    p.addlineofchange(descriptionsToAdd, prevToAdd, addcolor3);
    p.sequence();
    p.setoutputdir(tree.get<std::string>("general.outputdir"));
    p.subdetectors=tree.get<std::string>("general.subdetector");
    if(p.subdetectors=="")
    {
        p.subdetectors="BPIX,FPIX,TIB,TID,TOB,TEC,BPIX_Y,FPIX_Y";
    }
    p.drawparameter=tree.get<std::string>("general.parameter");
    if(p.drawparameter!="")
    {
        p.draw(p.subdetectors, p.drawparameter);
    }
    else
    {
        p.draw(p.subdetectors);
    }
}
