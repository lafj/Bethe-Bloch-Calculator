/*
*	bethe.cpp
*	
*	Created by Lars Fredrik Fjæra
*	Spring 2015
*	PHYS291 - Data Handling in Physics
*	UiB
*
*	This program is a graphical user interface "calculator" which
*	calculates and plots the energy loss and range of heavy charged
*	particles going through matter.
*	The program reads two files containing material- and particle
*	parameters. The user can choose what kind of particle and 
*	what material the particle should traverse. The user then defines
*	the energy of the incoming particle. Press 'Calculate and Draw' and the
*	program calculates the energy loss, the range and draws the plots.
*
*
*/

#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGComboBox.h>
#include <TGLabel.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include <TGTextEntry.h>
#include "iostream"
#include "Riostream.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TStyle.h"
#include "TObject.h"
#include "TTree.h"
#include "TMath.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TAxis.h"
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include "TGNumberEntry.h"
using namespace std;
void readfile();

class MyMainFrame {
    RQ_OBJECT("MyMainFrame")
private:
    TGMainFrame *fMain;
    TRootEmbeddedCanvas *fEcanvas;
    TGNumberEntry *userEnergy;
    TGNumberEntry *userCharge;
    TGNumberEntry *userMass;
    TGNumberEntry *userZ;
    TGNumberEntry *userA;
    TGNumberEntry *userIP;
    TGNumberEntry *userDens;
    TGComboBox *matList;
    TGComboBox *partList;
    TGLabel *answerlabel[4];
    TGLabel *calc[4];
    TGLabel *materiallabel[5];
    TGLabel *matparameterlabel[5];
    TGGroupFrame *resultframe;
    TGGroupFrame *userdefinedpart;
    TGGroupFrame *userdefinedmat;
    TGGroupFrame *materialframe;
    TGVerticalFrame *calcframe;
    TGVerticalFrame *matparameterframe;
public:
    MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
    virtual ~MyMainFrame();
    void DrawCalc();
};

struct Medium {
    string name;
    Double_t zoa;
    Double_t iPot;
    Double_t dens;
    Double_t atNo;
};
Medium *material[127];

struct IncomingPart {
    string partName;
    Int_t charge;
    Double_t mass;
};
IncomingPart *iPart[8];

struct Results {
    Double_t resEnLoss;
    Double_t resRange;
    Double_t resEnLossRho;
    Double_t resRangeRho;
};
Results output;

MyMainFrame::MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h) {
    // Create a main frame
    fMain = new TGMainFrame(p,w,h, kHorizontalFrame);
	
	// Create canvas widget
    fEcanvas = new TRootEmbeddedCanvas("Ecanvas",fMain,800,600);
    
    // Create a vertical frame widget for buttonframe and parameter inputs etc.
    TGVerticalFrame *leftframe = new TGVerticalFrame(fMain,200,40);
    
    //Label and ComboBox for incoming particle
    TGLabel *incPart = new TGLabel(leftframe, "\nChoose incoming particle:");
    leftframe->AddFrame(incPart, new TGLayoutHints(kLHintsTop));
    partList = new TGComboBox(leftframe);
    for (Int_t i=0;i<8;i++) {
        partList->AddEntry(iPart[i]->partName.c_str(),i);
    }
    partList->Select(0);
    partList->Resize(300,20);
    leftframe->AddFrame(partList, new TGLayoutHints(kLHintsTop));
    
    //Make a groupframe for user defined particle
    userdefinedpart = new TGGroupFrame(leftframe, "Define particle parametes",kVerticalFrame);
    leftframe->AddFrame(userdefinedpart, new TGLayoutHints(kLHintsLeft, 5,5,10,5));
    
    //Labels and input for user defined particle
    TGLabel *userchargelayout = new TGLabel(userdefinedpart, "\nCharge of incoming particle [e]:");
    userdefinedpart->AddFrame(userchargelayout);
    userCharge = new TGNumberEntry(userdefinedpart, 1, 9, 999, TGNumberFormat::kNESInteger,
    		TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 1, 1000);
    userCharge->Resize(100,20);
    userdefinedpart->AddFrame(userCharge);
    TGLabel *usermasslayout = new TGLabel(userdefinedpart, "\nMass of incoming particle [MeV/c^2]\n(Minimum 100 MeV/c^2):");
    userdefinedpart->AddFrame(usermasslayout);
    userMass = new TGNumberEntry(userdefinedpart, 100, 9, 999, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMin,100);
    userMass->Resize(100,20);
    userdefinedpart->AddFrame(userMass);
    
    //Label and ComboBox for material
    TGLabel *mater = new TGLabel(leftframe, "\nChoose material:");
    leftframe->AddFrame(mater, new TGLayoutHints(kLHintsTop));
    matList = new TGComboBox(leftframe);
    for (Int_t i=0;i<127;i++) {
        matList->AddEntry(material[i]->name.c_str(),i);
            }
    matList->Select(0);
    matList->Resize(300,20);
    leftframe->AddFrame(matList, new TGLayoutHints(kLHintsTop));

    //Make a frame for user defined material
	userdefinedmat = new TGGroupFrame(leftframe, "Define material parameters", kVerticalFrame);
	leftframe->AddFrame(userdefinedmat, new TGLayoutHints(kLHintsLeft, 5,5,10,5));
	TGLabel *userZlayout = new TGLabel(userdefinedmat, "\nAtomic nubmer (Z):");
	userdefinedmat->AddFrame(userZlayout);
	userZ = new TGNumberEntry(userdefinedmat,1,9,999,TGNumberFormat::kNESRealOne, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELNoLimits);
	userZ->Resize(100,20);
	userdefinedmat->AddFrame(userZ);
	TGLabel *userAlayout = new TGLabel(userdefinedmat, "\nMass number (A):");
	userdefinedmat->AddFrame(userAlayout);
	userA = new TGNumberEntry(userdefinedmat,1,9,999,TGNumberFormat::kNESRealOne, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELNoLimits);
	userA->Resize(100,20);
	userdefinedmat->AddFrame(userA);
	TGLabel *userIPlayout = new TGLabel(userdefinedmat, "\nIonisation Potential [eV]:");
	userdefinedmat->AddFrame(userIPlayout);
	userIP = new TGNumberEntry(userdefinedmat,1,9,999,TGNumberFormat::kNESRealOne, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELNoLimits);
	userIP->Resize(100,20);
	userdefinedmat->AddFrame(userIP);
	TGLabel *userDensLayout = new TGLabel(userdefinedmat, "\nDensity [g/cm^3]:");
	userdefinedmat->AddFrame(userDensLayout);
	userDens = new TGNumberEntry(userdefinedmat,1,9,999,TGNumberFormat::kNESRealOne, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELNoLimits);
	userDens->Resize(100,20);
	userdefinedmat->AddFrame(userDens);

    //Label and input for kinetic energy of incoming particle
    TGLabel *kinEnLabel = new TGLabel(leftframe, "\nKinetic energy of incoming particle [MeV]\n(From 1 to 10000 MeV):");
    leftframe->AddFrame(kinEnLabel, new TGLayoutHints(kLHintsTop));
    userEnergy = new TGNumberEntry(leftframe, 1, 9, 999, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 1, 10000);
    userEnergy->Resize(100,20);
    leftframe->AddFrame(userEnergy, new TGLayoutHints(kLHintsTop));
    
    // Create a horizontal frame for the buttons with Calc and Exit buttons
    TGHorizontalFrame *buttonframe = new TGHorizontalFrame(leftframe,200,40);
    leftframe->AddFrame(buttonframe, new TGLayoutHints(kLHintsTop | kLHintsCenterX));
    //Make Calc and Draw button
    TGTextButton *draw = new TGTextButton(buttonframe,"&Calcuate and draw");
    draw->Connect("Clicked()","MyMainFrame",this,"DrawCalc()");
    buttonframe->AddFrame(draw, new TGLayoutHints(kLHintsCenterX, 5,5,20,4));
    
    //Make Exit button
    TGTextButton *exit = new TGTextButton(buttonframe,"&Exit", "gApplication->Terminate(0)");
    buttonframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5,5,20,4));
    
    //Frame for outputting material parameters
    materialframe = new TGGroupFrame(leftframe,"Material parameters",kHorizontalFrame);
    materialframe->SetTitlePos(TGGroupFrame::kCenter);
    leftframe->AddFrame(materialframe, new TGLayoutHints(kLHintsTop,20,20,20,20));

    // Create a vertical frame widget that will display the results of the calculation
    resultframe = new TGGroupFrame(leftframe,"Results",kHorizontalFrame);
    resultframe->SetTitlePos(TGGroupFrame::kCenter);
    leftframe->AddFrame(resultframe, new TGLayoutHints(kLHintsBottom,20,20,1,20));
    
    TGVerticalFrame *matlabelframe = new TGVerticalFrame(materialframe,10,10);
    materialframe->AddFrame(matlabelframe, new TGLayoutHints(kLHintsLeft));

    matparameterframe = new TGVerticalFrame(materialframe,10,10);
    materialframe->AddFrame(matparameterframe,new TGLayoutHints(kLHintsRight));

    //Material parameter labels
    materiallabel[0] = new TGLabel(matlabelframe, "Type: ");
    matlabelframe->AddFrame(materiallabel[0]);
    materiallabel[1] = new TGLabel(matlabelframe, "Z: ");
    matlabelframe->AddFrame(materiallabel[1]);
    materiallabel[2] = new TGLabel(matlabelframe, "A: ");
    matlabelframe->AddFrame(materiallabel[2]);
    materiallabel[3] = new TGLabel(matlabelframe, "Ionisation Pot [eV]: ");
    matlabelframe->AddFrame(materiallabel[3]);
    materiallabel[4] = new TGLabel(matlabelframe, "Density [g/cm^3]: ");
    matlabelframe->AddFrame(materiallabel[4]);
    matparameterlabel[0] = new TGLabel(matparameterframe, "                     ");
    matparameterframe->AddFrame(matparameterlabel[0]);
    matparameterlabel[1] = new TGLabel(matparameterframe, "                     ");
    matparameterframe->AddFrame(matparameterlabel[1]);
    matparameterlabel[2] = new TGLabel(matparameterframe, "                     ");
    matparameterframe->AddFrame(matparameterlabel[2]);
    matparameterlabel[3] = new TGLabel(matparameterframe, "                     ");
    matparameterframe->AddFrame(matparameterlabel[3]);
    matparameterlabel[4] = new TGLabel(matparameterframe, "                     ");
    matparameterframe->AddFrame(matparameterlabel[4]);

    TGVerticalFrame *answerframe = new TGVerticalFrame(resultframe,10,10);
    resultframe->AddFrame(answerframe, new TGLayoutHints(kLHintsLeft));
    
    calcframe = new TGVerticalFrame(resultframe,10,10);
    resultframe->AddFrame(calcframe,new TGLayoutHints(kLHintsRight));
    
    //Result labels
    answerlabel[0] = new TGLabel(answerframe, "Energy loss [MeV/cm] = ");
    answerframe->AddFrame(answerlabel[0]);
    answerlabel[1] = new TGLabel(answerframe, "Energy loss [MeV cm^2/g] = ");
    answerframe->AddFrame(answerlabel[1]);
    answerlabel[2] = new TGLabel(answerframe, "Range [cm] = ");
    answerframe->AddFrame(answerlabel[2]);
    answerlabel[3] = new TGLabel(answerframe, "Range [g/cm^2] = ");
    answerframe->AddFrame(answerlabel[3]);
    calc[0] = new TGLabel(calcframe, "N/A            ");
    calcframe->AddFrame(calc[0]);
    calc[1] = new TGLabel(calcframe, "N/A            ");
    calcframe->AddFrame(calc[1]);
    calc[2] = new TGLabel(calcframe, "N/A            ");
    calcframe->AddFrame(calc[2]);
    calc[3] = new TGLabel(calcframe, "N/A            ");
    calcframe->AddFrame(calc[3]);

    // Adding frames
    fMain->AddFrame(leftframe, new TGLayoutHints(kLHintsLeft, 2,2,2,2));
    fMain->AddFrame(fEcanvas, new TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsExpandY, 10,10,10,10));
    
    // Set a name to the main frame
    fMain->SetWindowName("Bethe-Bloch Calculator");
    // Map all subwindows of main frame
    fMain->MapSubwindows();
    // Initialize the layout algorithm
    fMain->Resize(fMain->GetDefaultSize());
    // Map main frame
    fMain->MapWindow();
}

void MyMainFrame::DrawCalc() {
    
    Double_t e_in; //Input energy
    Int_t mat; //Material entryID variable
    Int_t partID; //Particle entryID variable
    Float_t k = 0.307; //constant in Bethe-Bloch
    Int_t z; //charge of incoming particle
    Double_t z_a; //atom number
    Double_t rho; //density of material
    Double_t beta,b2; //beta = v/c
    Double_t m_e = 0.511; //electron mass
    Double_t m_i; //mass of incoming particle
    Double_t gam_fact,g2; //gamma = 1/sqrt(1+beta^2)
    Double_t eta,eta2; //eta=beta*gam_fact
    Double_t t_max; //maximum energy transfer in single collision (own formula)
    Double_t ion_pot; //ionisation potential
    Double_t shell; //shell correction
    Double_t dEdx = 0; //Energy loss per unit length
    Double_t e_kin; //energy of incoming particle
    Double_t zdivA; // Z/A value for absorbing material
    Double_t A; // Absorbing material mass number
    
    //Getting user input
    e_in = userEnergy->GetNumber();
    mat = matList->GetSelected();
    partID = partList->GetSelected();
    
    // If user define their own incoming particle
    if (partID == 0) {
        z = userCharge->GetIntNumber();
        m_i = userMass->GetNumber();
    }
    else {
        z = iPart[partID]->charge;
        m_i = iPart[partID]->mass; //mEv
    }
	
    //If user define their own absorbing material
    if (mat == 0) {
        z_a = userZ->GetNumber();
	A = userA->GetNumber();
	rho = userDens->GetNumber();
	ion_pot = userIP->GetNumber()/1000000;
	zdivA = z_a/A;
    }
    else {
        z_a = material[mat]->atNo;
        zdivA = material[mat]->zoa;
        ion_pot = material[mat]->iPot/1000000;
        rho = material[mat]->dens;// g/cm^3
	A = pow(zdivA/z_a,-1);
    }


    //Output material parameters in window
    matparameterlabel[0]->SetText(material[mat]->name.c_str());
    matparameterlabel[1]->SetText(Form("%.1f",z_a));
    matparameterlabel[2]->SetText(Form("%.1f",A));
    matparameterlabel[3]->SetText(Form("%.1f",ion_pot*1000000));
    matparameterlabel[4]->SetText(Form("%E",rho));
   
    //Define arrays used for temporary calculations
    Double_t *eArray = new Double_t[10000];
    Double_t *kinArray = new Double_t[10000];
    Double_t *rArray = new Double_t[10000];
    
    //Defining different variables
    Int_t i = 0;
    Int_t nPoints1 = 0;
    Int_t nPoints = 0;
    Double_t rMax = 0;
    Double_t eMax = 0;
    Double_t xrMax = 0;
    Double_t step,dx;
    Double_t x = 0.000;
    Double_t dE;
    
    //Defining step length depending on energy
    dE = e_in/1000; //step size
    dx=0;
    e_kin=e_in;
    
    //Calculating output values
    while(e_kin>0 && dEdx>=0) {
        x = x+dx;
        
        gam_fact = (e_kin/m_i)+1;
        g2 = pow(gam_fact,2);
        b2 = 1-(1/g2);
        beta = sqrt(b2);
        eta = beta*gam_fact;
        eta2= g2*b2;
	
        shell = ((0.422377*pow(eta,-2)+0.0304043*pow(eta,-4)-0.00038106*pow(eta,-6))*pow(10,-6)*pow(ion_pot,2))+((3.850190*pow(eta,-2)+0.1667989*pow(eta,-4)-0.00157955*pow(eta,-6))*pow(10,-9)*pow(ion_pot,3));
        t_max = (2*m_e*eta2)/(1+(2*m_e/m_i)*sqrt(1+eta2)+(pow(m_e/m_i,2)));
        dEdx = k*pow(z,2)*zdivA*(1/b2)*(0.5*log(2*m_e*eta2*(t_max/pow(ion_pot,2)))-b2-shell/z_a);

        eArray[nPoints1] = dEdx;
        rArray[nPoints1] = x;
        
        dx = dE/dEdx;
        e_kin=e_kin-dE;

        nPoints1++;
    }
    
    //Outputting calculated values
    output.resEnLoss = eArray[0]*rho;
    output.resEnLossRho = eArray[0];
    //nPoits-1 because for nPoints either e_kin>0 or dEdx>=0 (or both)
    output.resRange = rArray[nPoints1-1]/rho;
    output.resRangeRho = rArray[nPoints1-1];

    //Defining arrays used for rangeplot
    Double_t *rPlot = new Double_t[10000];
    Double_t *xrPlot = new Double_t[10000];

	
	// New loop used to make rangeplot. Needed because stepsize is defined by path length of particle to make each
    // each point in the plot equally spaced.
	dEdx=0;
	dx = rArray[nPoints1-1]/1000;
	e_kin=e_in;
	x=0;
	Int_t nPoints2 = 0;
	
	while(e_kin>0 && dEdx>=0) {
        x = x+dx;
        
        gam_fact = (e_kin/m_i)+1;
        g2 = pow(gam_fact,2);
        b2 = 1-(1/g2);
        beta = sqrt(b2);
        eta = beta*gam_fact;
        eta2= g2*b2;
	
        shell = ((0.422377*pow(eta,-2)+0.0304043*pow(eta,-4)-0.00038106*pow(eta,-6))*pow(10,-6)*pow(ion_pot,2))+((3.850190*pow(eta,-2)+0.1667989*pow(eta,-4)-0.00157955*pow(eta,-6))*pow(10,-9)*pow(ion_pot,3));
        t_max = (2*m_e*eta2)/(1+(2*m_e/m_i)*sqrt(1+eta2)+(pow(m_e/m_i,2)));
        dEdx = k*pow(z,2)*zdivA*(1/b2)*(0.5*log(2*m_e*eta2*(t_max/pow(ion_pot,2)))-b2-shell/z_a);
        
        xrPlot[nPoints2] = dEdx;
        rPlot[nPoints2] = x;
                
        dE = dEdx*dx;
        e_kin=e_kin-dE;
	    
	    
        //If the increase between the two last energy loss calculations is more that ~7% of the previous value, the
        //loop stopps. This happens when the kinetic energy of the incoming particle apporaches zero.
        //This to make the bragg peak more defined. Doesn't affect range calculations.
        if (nPoints2>2) {
        	if (xrPlot[nPoints2]-xrPlot[nPoints2-1]>xrPlot[nPoints2-1]/15) {break;}
        }
        nPoints2++;
    }

	//Remove values where bethe bloch formula fails.
	rPlot[nPoints2] = rPlot[nPoints2-1];
    	xrPlot[nPoints2] = 0;
   	xrPlot[nPoints2-1] = 0;
	
   	//Finding maximumvalues used to define axis.
	for (i=0;i<nPoints2+1;i++) {
		if (rPlot[i]>rMax) {
            	rMax = rPlot[i];
        	}
        	if (xrPlot[i]>xrMax) {
            	xrMax=xrPlot[i];
        	}
    }
    
    
    //Yet another loop of calculatios for energy loss plot. Incoming energy for the plot is always the same.
	//It makes the plot look nicer.
    
	e_kin = 10000.0; // define e_kin for plotting dE/dx
    step = e_kin/1000.0000;
    nPoints=0;
    
    for(i=0;i<1000;i++) {
        gam_fact = (e_kin/m_i)+1;
        g2 = pow(gam_fact,2);
        b2 = 1-(1/g2);
        beta = sqrt(b2);
        eta = beta*gam_fact;
        eta2= g2*b2;
        //Stop the loop for low energy
        if (eta <= 0.1) {break;}
        shell = ((0.422377*pow(eta,-2)+0.0304043*pow(eta,-4)-0.00038106*pow(eta,-6))*pow(10,-6)*pow(ion_pot,2))+((3.850190*pow(eta,-2)+0.1667989*pow(eta,-4)-0.00157955*pow(eta,-6))*pow(10,-9)*pow(ion_pot,3));
        t_max = (2*m_e*eta2)/(1+(2*m_e/m_i)*sqrt(1+eta2)+(pow(m_e/m_i,2)));
        dEdx = k*pow(z,2)*zdivA*(1/b2)*(0.5*log(2*m_e*eta2*(t_max/pow(ion_pot,2)))-b2-shell/z_a);
        
        eArray[i] = dEdx;
        kinArray[i] = e_kin;
        
        e_kin=e_kin-step;
        nPoints++;
    }
    //Arrays for energy plot
    Double_t *ePlot = new Double_t[nPoints];
    Double_t *dedxPlot = new Double_t[nPoints];
    
    //Fill plotting arrays
    for (i=0;i<nPoints;i++) {
        ePlot[i] = kinArray[i];
        dedxPlot[i] = eArray[i];
        //Find max value
        if (dedxPlot[i]>eMax) {
            eMax = dedxPlot[i];
        }
    }
    //Canvas and pads
    TCanvas *c1 = fEcanvas->GetCanvas();
    c1->cd();
    c1->SetFillColor(10);
    c1->Modified();
    c1->Update();
    TPad *pad1 = new TPad("pad1","Energy Loss",0.02,0.51,0.98,0.98,10);
    TPad *pad2 = new TPad("pad2","Range",0.02,0.02,0.98,0.49,10);
    pad1->Draw();
    pad2->Draw();
    
    //Energy loss plot
    pad1->cd();
    TGraph *gr1 = new TGraph(nPoints,ePlot,dedxPlot);
    pad1->SetLogy();
    pad1->SetLogx();
    gr1->GetYaxis()->SetRangeUser(0,2*eMax);
    gr1->GetXaxis()->SetLimits(0,e_kin);
    gr1->GetXaxis()->SetTitle("\nKinetic energy [MeV]");
    gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->SetTitle("dE/dx [MeV cm^{2}/g]");
    gr1->GetYaxis()->CenterTitle();
    gr1->SetTitle("Energy Loss");
    gr1->Draw("APC");
    c1->Update();
    
    //Range plot
    pad2->cd();
    TGraph *gr2 = new TGraph(nPoints2+1,rPlot,xrPlot);
    gr2->GetYaxis()->SetRangeUser(0,xrMax+xrMax/5);
    gr2->GetXaxis()->SetLimits(0,rMax+rMax/10);
    gr2->GetXaxis()->SetTitle("Range [g/cm^{2}]");
    gr2->GetXaxis()->CenterTitle();
    gr2->GetYaxis()->SetTitle("dE/dx [MeV cm^{2}/g]");
    gr2->GetYaxis()->CenterTitle();
    gr2->SetTitle("Particle Range");
    gr2->Draw("APC");
    c1->Update();
    
    //Output results in window
    calc[0]->SetText(Form("%E",output.resEnLoss));
    calc[1]->SetText(Form("%E",output.resEnLossRho));
    calc[2]->SetText(Form("%E",output.resRange));
    calc[3]->SetText(Form("%E",output.resRangeRho));
    calcframe->Layout();
}

void readfile() {
	// initialize the structs
	for( Int_t zxc=0;zxc< 127;zxc++) material[zxc]=new Medium;
	for( Int_t zxc=0;zxc< 8;zxc++)  iPart[zxc] = new IncomingPart;

    // Read material file
    string line;
    Int_t count = 0;
    string totArray[127*5];
    ifstream in;
    in.open("materials.txt");
    while (getline(in,line)) {
        istringstream iss(line);
        string field;
        while (getline(iss, field, '\t')){
            totArray[count] = field;
            count++;
        }
    }
    //Filling the structs
    Int_t i;
    count=0;
    for (i=0;i<127*5;i=i+5) {
        material[count]->name = totArray[i];
        count++;
    }
    count=0;
    for (i=1;i<127*5;i=i+5) {
        material[count]->zoa = atof(totArray[i].c_str() );
        count++;
    }
    count=0;
    for (i=2;i<127*5;i=i+5) {
        material[count]->iPot = atof(totArray[i].c_str() );
        count++;
    }
    count=0;
    for (i=3;i<127*5;i=i+5) {
        material[count]->dens = atof(totArray[i].c_str() );
        count++;
    }
    count=0;
    for (i=4;i<127*5;i=i+5) {
        material[count]->atNo = atof(totArray[i].c_str() );
        count++;
    }
    in.close();
    
    //Read particle file
    count = 0;
    string line1;
    string partArray[8*3];
    ifstream in1;
    in1.open("particles.txt");
    while (getline(in1,line1)) {
        istringstream iss(line1);
        string field1;
        while (getline(iss, field1, '\t')){
            partArray[count] = field1;
            count++;
        }
    }
    //Filling the structs
    count=0;
    for (i=0;i<8*3;i=i+3) {
        iPart[count]->partName = partArray[i];
        count++;
        }
    count=0;
    for (i=1;i<8*3;i=i+3) {
        iPart[count]->charge = atof(partArray[i].c_str());
        count++;
        }
    count=0;
    for (i=2;i<8*3;i=i+3) {
        iPart[count]->mass = atof(partArray[i].c_str());
        count++;
        }
    in1.close();
}

MyMainFrame::~MyMainFrame() {
    // Clean up used widgets: frames, buttons, layout hints
    fMain->Cleanup();
    delete fMain;}

void bethe() { //Main function
    readfile();
    new MyMainFrame(gClient->GetRoot(),1000,1000);
}


