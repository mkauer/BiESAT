{
gROOT->Reset();
gROOT->SetStyle("Plain");

c1 = new TCanvas("c1","Results",700,500);
c1->SetFillColor(0);
c1->SetGrid();
c1->GetFrame()->SetFillColor(0);
c1->GetFrame()->SetBorderSize(12);

TH1F *mg = c1->DrawFrame(4,6,22,15);
mg->SetXTitle("Scintillator Size (cm)");
mg->SetYTitle("Resolution FWHM @ 1MeV");
mg->SetTitle("8in PMT's with Lightguide Results");


Color_t color1=8;
Color_t color2=46;
Color_t color3=9;

Int_t n1=2;
Float_t x1[n1]={5,14};
Float_t y1[n1]={8.5,9.2};
Float_t ex1[n1]={};
Float_t ey1[n1]={8.5*.03,9.2*.03};
etl = new TGraphErrors(n1,x1,y1,ex1,ey1);
etl->SetMarkerColor(color1);
etl->SetMarkerStyle(24);
etl->SetMarkerSize(2);
etl->Draw("P");

etltex1 = new TText(5.5,9,"Teflon");
etltex1->SetTextColor(color1);
etltex1->SetTextSize(.05);
etltex1->SetTextAngle(45);
etltex1->Draw();

etltex2 = new TText(14.5,8.8,"Teflon");
etltex2->SetTextColor(color1);
etltex2->SetTextSize(.05);
etltex2->SetTextAngle(45);
etltex2->Draw();


Int_t n2=1;
Float_t x2[n2]={5};
Float_t y2[n2]={7.9};
Float_t ex2[n2]={};
Float_t ey2[n2]={7.9*.03};
ham = new TGraphErrors(n2,x2,y2,ex2,ey2);
ham->SetMarkerColor(color2);
ham->SetMarkerStyle(25);
ham->SetMarkerSize(2);
ham->Draw("P");

hamtex1 = new TText(5.5,7.7,"Teflon");
hamtex1->SetTextColor(color2);
hamtex1->SetTextSize(.05);
hamtex1->SetTextAngle(45);
hamtex1->Draw();


Int_t n3=7;
Float_t x3[n3]={14,15,9,5,14,5,20};
Float_t y3[n3]={11.5,10.4,10.1,8.5,12.7,7.5,10.9};
Float_t ex3[n3]={};
Float_t ey3[n3]={11.5*.03,10.4*.03,10.1*.03,8.5*.03,12.7*.03,7.5*.03,10.9*.03};
sba = new TGraphErrors(n3,x3,y3,ex3,ey3);
sba->SetMarkerColor(color3);
sba->SetMarkerStyle(26);
sba->SetMarkerSize(2);
sba->Draw("P");

sbatex1 = new TText(5.5,7,"Teflon");
sbatex1->SetTextColor(color3);
sbatex1->SetTextSize(.05);
sbatex1->SetTextAngle(45);
sbatex1->Draw();

sbatex2 = new TText(5.5,8.3,"R-ESR");
sbatex2->SetTextColor(color3);
sbatex2->SetTextSize(.05);
sbatex2->SetTextAngle(45);
sbatex2->Draw();

sbatex3 = new TText(14.5,12.5,"Teflon");
sbatex3->SetTextColor(color3);
sbatex3->SetTextSize(.05);
sbatex3->SetTextAngle(45);
sbatex3->Draw();

sbatex4 = new TText(14.5,11.5,"Mylar");
sbatex4->SetTextColor(color3);
sbatex4->SetTextSize(.05);
sbatex4->SetTextAngle(45);
sbatex4->Draw();

sbatex5 = new TText(15.5,10.5,"Mylar");
sbatex5->SetTextColor(color3);
sbatex5->SetTextSize(.05);
sbatex5->SetTextAngle(45);
sbatex5->Draw();

sbatex6 = new TText(9.5,9.8,"Mylar");
sbatex6->SetTextColor(color3);
sbatex6->SetTextSize(.05);
sbatex6->SetTextAngle(45);
sbatex6->Draw();

sbatex7 = new TText(20.5,11.1,"R-ESR");
sbatex7->SetTextColor(color3);
sbatex7->SetTextSize(.05);
sbatex7->SetTextAngle(45);
sbatex7->Draw();


TLegend *leg = new TLegend(0.7,0.85,0.995,0.995);
leg->SetFillColor(0);
leg->AddEntry(etl,"ETL PM","p");
leg->AddEntry(ham,"Hamamatsu PM","p");
leg->AddEntry(sba,"Ham. SBA PM","p");
leg->Draw();

c1->Print("8in_PMs.png");
}

