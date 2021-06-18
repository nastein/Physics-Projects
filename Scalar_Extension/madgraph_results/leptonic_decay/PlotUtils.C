
//============================================================

class Plot {

public:
    Plot(string id);
    void Add(string leg, int color, int style, TH1D* ph, bool scale = false, double scalefactor = 0);
    void Add(string leg, int color, int style, TProfile* ph);
    void Add(string leg, int color, int style, int fillstyle, TH1D* ph);
    void Add(string leg, int color, int style, int fillstyle, TH2D* ph);
    void Add(string leg, int color, TGraph* ph);
    void Add(string leg, int color, TGraphErrors* ph, vector<Double_t> y_errs);
    void Title(string t) {title_ = t;};
    void Draw(string varX, string varY, const double posX = 0.16, const double scale_axis = 1.1, const double min_y = 0);
    void DrawDiff(string varX, string varY, const double posX = 0.16, const double scale_axis = 1.1, const double min_y = 0);
    void DrawRatio(string varX, string varY, const double posX = 0.16, const double scale_axis = 1.1, const double min_y = 0);
    void DrawGraphRatio(string varX, string varY, const double posX = 0.16, int colors[] = NULL);
    void Draw2D(string varX, string varY, const double posX = 0.16, const double scale_axis = 1.1);
    void DrawProfile(string varX, string varY, const double posX, const double scale_axis);
    void DrawGraph(string varX, string varY, const double posX);
    void SetRedHeatPalette() const;
    void AddPlotLabel(const char* label, const double x, const double y, const double size = 0.05);
    void Drawline(double xmin, double ymin, double xmax, double ymax);
    void SetLog(bool b=true) {log_ = b;};

    double ratio_err(double x, double y, double x_err, double y_err) {
        return (1/y)*sqrt(x_err*x_err + pow(x/y,2)*y_err*y_err);
    }

private:
    string id_;
    int num_items_;
    vector<TH1D*> Analysis_hist_;
    vector<string> Analysis_leg_;
    vector<TH2D*> Analysis_hist2D_;
    vector<TProfile*> Analysis_TPr_;
    vector<TGraph*> Analysis_graph_;
    bool log_;
    string title_;
};


Plot::Plot(string id) {
    id_ = id;
    title_ = "";
    num_items_ = 0;
    Analysis_hist_.clear();
    Analysis_hist2D_.clear();
    Analysis_TPr_.clear();
    Analysis_graph_.clear();
    Analysis_leg_.clear();

    log_ = false;
}


void Plot::Add(string leg, int color, int style, TH1D* ph, bool scale, double scalefactor) {

    num_items_++;
    TH1D* h = ph;
	h->SetLineColor(color);
	h->SetLineStyle(style);
    if (scale == true) {
        h->Sumw2();
        h->Scale(scalefactor);
    }
	Analysis_hist_.push_back(h);
	Analysis_leg_.push_back(leg);

}


void Plot::Add(string leg, int color, int style, TProfile* ph) {

    num_items_++;
    TProfile* h = ph;
    h->SetLineColor(color);
    h->SetLineStyle(style);
    Analysis_TPr_.push_back(h);
    Analysis_leg_.push_back(leg);

}


void Plot::Add(string leg, int color, int style, int fillystle, TH2D* ph) {

    num_items_++;
    TH2D* h = ph;
    h->SetFillColor(color);
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    Analysis_hist2D_.push_back(h);
    Analysis_leg_.push_back(leg);


}

void Plot::Add(string leg, int color, TGraph* ph) {

    num_items_++;
    TGraph *g = ph;
    g->SetLineColor(color);
    g->SetLineWidth(3);
    //g->SetMarkerColor(color);
    //g->SetMarkerStyle(21);
    Analysis_graph_.push_back(g);
    Analysis_leg_.push_back(leg);


}

void Plot::Draw(string varX, string varY, const double posX, const double scale_axis, const double min_y) {

    if (num_items_ == 0) return;

    TCanvas *can1 = new TCanvas((""+id_).c_str(), (""+id_).c_str(), 700, 500);
    if (log_) can1->SetLogy();

    TLegend *leg = new TLegend(posX, .7, posX+0.20, 0.9, "L NDC");
    leg->SetTextAlign(11);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(45);//was 45
    leg->SetTextSize(20);

    double maxY = 0.0;
    for (unsigned int i = 0; i < num_items_; i++) {
    	TH1D* h = Analysis_hist_[i];
    	if (h->GetMaximum() > maxY) {maxY = h->GetMaximum();}
    }

    for (unsigned int i = 0; i < num_items_; i++) {
    	TH1D* h = (TH1D*) Analysis_hist_[i]->Clone();

        h->SetLineWidth(2);

        if (log_) h->SetMinimum(0.001);
        else if (min_y !=0 ) h->SetMinimum(min_y);
        else h->SetMinimum(0);

    	if (i == 0) {

    	    h->SetMaximum(scale_axis*maxY);  /// \todo : this is not well done, don't change h
            if (varY != "") h->SetYTitle(varY.c_str());
    	    h->SetXTitle(varX.c_str());
    	    h->GetXaxis()->CenterTitle(-1);
    	    h->GetYaxis()->CenterTitle(-1);
            h->GetXaxis()->SetLabelFont(42);
            h->GetYaxis()->SetLabelFont(42);
    	    h->SetNdivisions(504, "y");
    	    h->SetNdivisions(504, "x");
    	    h->SetTitle(title_.c_str());
    	    gStyle->SetTitleFontSize(0.07);
    	    h->Draw("hist E");
    	} else {
            h->SetMaximum(scale_axis*maxY);  /// \todo : this is not well done, don't change h
            if (varY != "") h->SetYTitle(varY.c_str());
            h->SetXTitle(varX.c_str());
            h->GetXaxis()->CenterTitle(-1);
            h->GetYaxis()->CenterTitle(-1);
            h->GetXaxis()->SetLabelFont(42);
            h->GetYaxis()->SetLabelFont(42);
            h->SetNdivisions(504, "y");
            h->SetNdivisions(504, "x");
            h->SetTitle(title_.c_str());
            gStyle->SetTitleFontSize(0.07);
            h->Draw("hist same E");
    	}
	   if (!Analysis_leg_[i].empty()) {leg->AddEntry(h, Analysis_leg_[i].c_str());};
    }

    leg->Draw();

}

void Plot::DrawDiff(string varX, string varY, const double posX, const double scale_axis, const double min_y) {

    //The first histogram added must be the base histogram

    if (num_items_ < 2) return;

    TCanvas *can1 = new TCanvas((""+id_).c_str(), (""+id_).c_str(), 700, 500);
    if (log_) can1->SetLogy();

    TLegend *leg = new TLegend(posX, .7, posX+0.20, 0.9, "L NDC");
    leg->SetTextAlign(11);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(45);//was 45
    leg->SetTextSize(20);

    /*
    double maxY = 0.0;
    for (unsigned int i = 0; i < num_items_; i++) {
        TH1D* h = Analysis_hist_[i];
        if (h->GetMaximum() > maxY) {maxY = h->GetMaximum();}
    }
    */

    TH1D* base_h = (TH1D*) Analysis_hist_[0]->Clone();

    for (unsigned int i = 1; i < num_items_; i++) {
        TH1D* h = (TH1D*) Analysis_hist_[i]->Clone();
        h->Add(base_h, -1.);

        h->SetLineWidth(2);

        if (log_) h->SetMinimum(0.001);
        else if (min_y !=0 ) h->SetMinimum(min_y);
        else h->SetMinimum(0);

        if (i == 1) {

            //h->SetMaximum(scale_axis*maxY);  /// \todo : this is not well done, don't change h
            if (varY != "") h->SetYTitle(varY.c_str());
            h->SetXTitle(varX.c_str());
            h->GetXaxis()->CenterTitle(-1);
            h->GetYaxis()->CenterTitle(-1);
            h->GetXaxis()->SetLabelFont(42);
            h->GetYaxis()->SetLabelFont(42);
            h->SetNdivisions(504, "y");
            h->SetNdivisions(504, "x");
            h->SetTitle(title_.c_str());
            gStyle->SetTitleFontSize(0.07);
            h->Draw("hist E");
        } else {
            if (varY != "") h->SetYTitle(varY.c_str());
            h->SetXTitle(varX.c_str());
            h->GetXaxis()->CenterTitle(-1);
            h->GetYaxis()->CenterTitle(-1);
            h->GetXaxis()->SetLabelFont(42);
            h->GetYaxis()->SetLabelFont(42);
            h->SetNdivisions(504, "y");
            h->SetNdivisions(504, "x");
            h->SetTitle(title_.c_str());
            gStyle->SetTitleFontSize(0.07);
            h->Draw("hist same E");
        }
       if (!Analysis_leg_[i].empty()) {leg->AddEntry(h, Analysis_leg_[i].c_str());};
    }

    leg->Draw();

}

void Plot::Draw2D(string varX, string varY, const double posX, const double scale_axis) {

    if (num_items_ == 0) return;

    TCanvas *can1 = new TCanvas((""+id_).c_str(), (""+id_).c_str(), 700, 500);
    if (log_) can1->SetLogy();
//    can1->SetRightMargin(0.15);

    TLegend *leg = new TLegend(posX, 0.98 - 0.06*(3+num_items_), posX+0.20, 0.98, "L NDC");
    leg->SetTextAlign(11);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(45);
    leg->SetTextSize(30);

    double maxY = 0.0;
    for (unsigned int i = 0; i < num_items_; i++) {
        TH2D* h = Analysis_hist2D_[i];
        if (h->GetMaximum() > maxY) {maxY = h->GetMaximum();}
    }

    for (unsigned int i = 0; i < num_items_; i++) {
        TH2D* h = (TH2D*) Analysis_hist2D_[i]->Clone();

    //  h->SetLineStyle(1);
        //h->SetLineWidth(2);

        if (i == 0) {
            h->SetMaximum(scale_axis*maxY);  /// \todo : this is not well done, don't change h
            h->SetYTitle(varY.c_str());
            h->SetXTitle(varX.c_str());
            h->GetXaxis()->CenterTitle(-1);
            h->GetYaxis()->CenterTitle(-1);
            h->SetNdivisions(504, "y");
            h->SetNdivisions(504, "x");
            h->SetTitle(title_.c_str());
            gStyle->SetTitleFontSize(0.03);
            h->Draw("COLZ");
        } else {
            h->Draw("SAME COLZ");
        }
       if (!Analysis_leg_[i].empty()) {leg->AddEntry(h, Analysis_leg_[i].c_str());};
    }

    leg->Draw();

    //can1->Print(("MUSE_2d_vertex_"+id_+".pdf").c_str());
}

void Plot::DrawProfile(string varX, string varY, const double posX, const double scale_axis) {

    if (num_items_ == 0) return;

    TCanvas *can1 = new TCanvas((""+id_).c_str(), (""+id_).c_str(), 700, 500);
    if (log_) can1->SetLogy();
//    can1->SetRightMargin(0.15);

    TLegend *leg = new TLegend(posX, 0.98 - 0.06*(3+num_items_), posX+0.20, 0.98, "L NDC");
    leg->SetTextAlign(11);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(45);
    leg->SetTextSize(30);

    double maxY = 0.0;
    for (unsigned int i = 0; i < num_items_; i++) {
        TProfile* h = Analysis_TPr_[i];
        if (h->GetMaximum() > maxY) {maxY = h->GetMaximum();}
    }

    for (unsigned int i = 0; i < num_items_; i++) {
        TProfile* h = (TProfile*) Analysis_TPr_[i]->Clone();

    //  h->SetLineStyle(1);
        //h->SetLineWidth(2);

        if (i == 0) {
            h->SetMaximum(scale_axis*maxY);  /// \todo : this is not well done, don't change h
            h->SetYTitle(varY.c_str());
            h->SetXTitle(varX.c_str());
            h->GetXaxis()->CenterTitle(-1);
            h->GetYaxis()->CenterTitle(-1);
            h->SetNdivisions(504, "y");
            h->SetNdivisions(504, "x");
            h->SetTitle(title_.c_str());
            gStyle->SetTitleFontSize(0.03);
            h->Draw("E hist");
        } else {
            h->Draw("E hist same");
        }
       if (!Analysis_leg_[i].empty()) {leg->AddEntry(h, Analysis_leg_[i].c_str());};
    }

    leg->Draw();

    //can1->Print(("MUSE_2d_vertex_"+id_+".pdf").c_str());
}

void Plot::DrawGraphRatio(string varX, string varY, const double posX, int colors[]) {

    if (num_items_ < 2) return;
    TCanvas *can1 = new TCanvas((id_).c_str(), (id_).c_str(), 800, 800);

    auto *p2 = new TPad("p2","p2",0.,0.,1.,0.3); p2->Draw();
    p2->SetTopMargin(0.001);
    p2->SetBottomMargin(0.3);
    p2->SetGrid();

    auto *p1 = new TPad("p1","p1",0.,0.3,1.,1.);  p1->Draw();
    p1->SetBottomMargin(0.001);
    p1->cd();
    p1->SetGrid();

    TLegend *leg = new TLegend(posX, 0.98 - 0.06*(3+num_items_), posX+0.20, 0.98, "L NDC");
    leg->SetTextAlign(11);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(45);//was 45
    leg->SetTextSize(20);

    double maxY = 0.0;
    for (unsigned int i = 0; i < num_items_; i++) {
        TGraph* h = Analysis_graph_[i];
        int n = h->GetN();
        double* y = h->GetY();
        int locmax = TMath::LocMax(n,y);
        double tmax = y[locmax];
        if (tmax > maxY) {maxY = tmax;}
    }

    TGaxis::SetMaxDigits(3);
    for (unsigned int i = 0; i < num_items_; i++) {
        TGraph* h = (TGraph*) Analysis_graph_[i]->Clone();


        if (i == 0) {

            h->GetYaxis()->SetRangeUser(0,maxY*1.2);  /// \todo : this is not well done, don't change h
            //h->SetMaximum(1.3);
            if (log_) h->SetMinimum(0.01);
            else h->SetMinimum(0.0);
            h->GetYaxis()->SetTitle(varY.c_str());
            h->GetXaxis()->SetTitle(varX.c_str());
            h->GetXaxis()->CenterTitle(-1);
            //h->SetMarkerStyle(21);
            //h->SetMark
            h->GetYaxis()->CenterTitle(-1);
            h->SetTitle(title_.c_str());
            gStyle->SetTitleFontSize(0.03);
            h->Draw();
        } else {
            h->Draw("same");
        }
       if (!Analysis_leg_[i].empty()) {leg->AddEntry(h, Analysis_leg_[i].c_str(), "l");};
    }

    leg->Draw();

    p2->cd();
    TGraph *r[num_items_ - 1];
    double x,y;//graph[0]
    double z,w;//graph[i]



    for (int i = 0; i < num_items_ - 1; i++) {
      r[i] = new TGraph(); r[i]->SetTitle(title_.c_str());
      for (int j = 0; j < Analysis_graph_[0]->GetN(); j++) {
        Analysis_graph_[0]->GetPoint(j,x,y);
        Analysis_graph_[i+1]->GetPoint(j,z,w);
        r[i]->SetPoint(j, x, w/y);
      }

      r[i]->GetXaxis()->SetLabelSize(0.075);
      r[i]->GetYaxis()->SetLabelSize(0.075);
      r[i]->SetLineColor(colors[i+1]);
      r[i]->SetLineWidth(3);
      //r[i]->SetMarkerColor(colors[i+1]);
      //r[i]->SetMarkerStyle(21);
    }

    r[0]->Draw("");

    for (int i = 1; i < num_items_ - 1; i++) {
        r[i]->Draw("same");
    }

}

void Plot::DrawRatio(string varX, string varY, const double posX, const double scale_axis, const double min_y) {


    if (num_items_ == 0) return;

    TCanvas *can1 = new TCanvas((""+id_).c_str(), (""+id_).c_str(), 800, 800);
    if (log_) can1->SetLogy();
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.45, 1, 1);
    pad1->SetBottomMargin(0);
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    TLegend *leg = new TLegend(posX, 0.98 - 0.06*(3+num_items_), posX+0.20, 0.98, "L NDC");
    leg->SetTextAlign(11);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(45);
    leg->SetTextSize(18);

    TH1D* h[num_items_];
    TH1D* h_ratio[num_items_];

    double maxY = 0.0;
    for (unsigned int i = 0; i < num_items_; i++) {
        TH1D* h = Analysis_hist_[i];
        if (h->GetMaximum() > maxY) {maxY = h->GetMaximum();}
    }

    for (unsigned int k =0; k <= 1; k++) {

        if (k==1) {
            can1->cd();
            TPad *pad2 = new TPad("pad2", "pad2", 0, 0.1, 1, 0.45);
            pad2->SetTopMargin(0);
            pad2->SetBottomMargin(.2);
            pad2->Draw();
            pad2->cd();
        }
        for (unsigned int i =0; i < num_items_; i++) {

            if (k == 0) {
                h[i] = Analysis_hist_[i];
                (h[i])->SetLineWidth(2);

                if (i == 0) {

                    h[i]->SetMaximum(scale_axis*maxY);  /// \todo : this is not well done, don't change h
                    if (log_) h[i]->SetMinimum(0.000001);
                    else if (min_y !=0 ) h[i]->SetMinimum(min_y);
                    else h[i]->SetMinimum(0.0);
                    h[i]->SetYTitle(varY.c_str());
                    h[i]->GetXaxis()->CenterTitle(-1);
                    h[i]->GetYaxis()->CenterTitle(-1);
                    h[i]->SetNdivisions(504, "y");
                    h[i]->SetNdivisions(504, "x");
                    h[i]->SetTitle(title_.c_str());
                    gStyle->SetTitleFontSize(0.03);
                    h[i]->Draw("hist E");
                } else {
                    h[i]->DrawCopy("hist E same");
                }
                if (!Analysis_leg_[i].empty()) {leg->AddEntry(h[i], Analysis_leg_[i].c_str());};

                leg->Draw();
            }

            else {

                if (i == 0) continue;

                h_ratio[i] = (TH1D*) h[i]->Clone();
                h_ratio[i]->SetMaximum(2.0);
                h_ratio[i]->Sumw2();
                h_ratio[i]->Divide(h[0]);
                h_ratio[i]->SetMinimum(0.8);
                h_ratio[i]->Draw("hist E same");
                h_ratio[i]->SetTitle("");
                h_ratio[i]->GetYaxis()->SetTitle("Ratio");
                h_ratio[i]->GetYaxis()->SetNdivisions(505);
                h_ratio[i]->GetYaxis()->SetTitleSize(20);
                h_ratio[i]->GetYaxis()->SetTitleFont(43);
                h_ratio[i]->GetYaxis()->SetTitleOffset(1.55);
                h_ratio[i]->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
                h_ratio[i]->GetYaxis()->SetLabelSize(15);

                h_ratio[i]->GetXaxis()->SetTitleSize(20);
                h_ratio[i]->GetXaxis()->SetTitleFont(43);
                h_ratio[i]->GetXaxis()->SetTitleOffset(4.);
                h_ratio[i]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
                h_ratio[i]->GetXaxis()->SetLabelSize(15);
                h_ratio[i]->GetXaxis()->SetTitle(varX.c_str());

            }
        }
    }

    pad1->cd();
}


void Plot::DrawGraph(string varX, string varY, const double posX) {

    if (num_items_ == 0) return;

    TCanvas *can1 = new TCanvas((""+id_).c_str(), (""+id_).c_str(), 700, 500);
    if (log_) can1->SetLogy();
//    can1->SetRightMargin(0.15);

    TLegend *leg = new TLegend(posX, 0.98 - 0.06*(3+num_items_), posX+0.20, 0.98, "L NDC");
    leg->SetTextAlign(11);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(45);//was 45
    leg->SetTextSize(20);

    double maxY = 0.0;
    for (unsigned int i = 0; i < num_items_; i++) {
        TGraph* h = Analysis_graph_[i];
        int n = h->GetN();
        double* y = h->GetY();
        int locmax = TMath::LocMax(n,y);
        double tmax = y[locmax];
        if (tmax > maxY) {maxY = tmax;}
    }

    TGaxis::SetMaxDigits(3);

    for (unsigned int i = 0; i < num_items_; i++) {
        TGraph* h = (TGraph*) Analysis_graph_[i]->Clone();


        if (i == 0) {

            h->GetYaxis()->SetRangeUser(0,maxY*1.2);  /// \todo : this is not well done, don't change h
            //h->SetMaximum(1.3);
            if (log_) h->SetMinimum(0.01);
            else h->SetMinimum(0.0);
            h->GetYaxis()->SetTitle(varY.c_str());
            h->GetXaxis()->SetTitle(varX.c_str());
            h->GetXaxis()->CenterTitle(-1);
            h->GetYaxis()->CenterTitle(-1);
            h->SetTitle(title_.c_str());
            gStyle->SetTitleFontSize(0.03);
            h->Draw();
        } else {
            h->Draw("same");
        }
       if (!Analysis_leg_[i].empty()) {leg->AddEntry(h, Analysis_leg_[i].c_str(), "l");};
    }

    leg->Draw();

    //can1->Print(("MUSE_Z_vertex_"+id_+".pdf").c_str());

}

void Plot::SetRedHeatPalette() const
{
     const int NRGBs = 9;
     static bool initialized=false;
     const int n_color_contours = 999;
     static int* colors=new int[n_color_contours];

    if(!initialized){
        // White -> red
        Double_t stops[NRGBs] = { 0.00, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000};
        Double_t red[NRGBs]   = { 1.00, 1.00, 0.99, 0.99, 0.98, 0.94, 0.80, 0.65, 0.40 };
        Double_t green[NRGBs] = { 0.96, 0.88, 0.73, 0.57, 0.42, 0.23, 0.09, 0.06, 0.00 };
        Double_t blue[NRGBs]  = { 0.94, 0.82, 0.63, 0.45, 0.29, 0.17, 0.11, 0.08, 0.05 };
        int colmin=TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, n_color_contours);
        for(uint i=0; i<n_color_contours; ++i) colors[i]=colmin+i;

        initialized=true;
     }
     gStyle->SetNumberContours(n_color_contours);
     gStyle->SetPalette(n_color_contours, colors);
}

void Plot::AddPlotLabel(const char* label, const double x, const double y, const double size)
{
    TLatex *latex = new TLatex( x, y, label );

    latex->SetNDC();
    latex->SetTextSize(size);
    latex->SetTextColor(1);
    latex->SetTextFont(52);
    latex->SetTextAlign(22);
    latex->SetTextAngle(0);
    latex->Draw();
}

void Plot::Drawline(double xmin, double ymin, double xmax, double ymax)
{
    TLine *line = new TLine(xmin, ymin, xmax, ymax);
    line->SetLineColor(kRed);
    //line->SetLineWidth(.5);
    line->Draw();
}

// ==========================================================