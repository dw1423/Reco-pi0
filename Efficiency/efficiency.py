import ROOT
ROOT.gStyle.SetOptStat(0)
from array import array

name_str = ['miniTreeAM_modifEcal2_low.root',
            'miniTreeAM_modifEcal1p5_low.root',
            'miniTreeAM_modifEcal1_low.root',
            'miniTreeAM_def_low.root'   
]
file = ROOT.TFile(name_str[0])
tree = file.Get("outtree")

# Define vectors for branches
pho_e = ROOT.std.vector('double')()
pho_px = ROOT.std.vector('double')()
pho_py = ROOT.std.vector('double')()
pho_pz = ROOT.std.vector('double')()
genpho_e = ROOT.std.vector('double')()
genpho_px = ROOT.std.vector('double')()
genpho_py = ROOT.std.vector('double')()
genpho_pz = ROOT.std.vector('double')()
genpi0_e = ROOT.std.vector('double')()
genpi0_m = ROOT.std.vector('double')()

# Set branch addresses
tree.SetBranchAddress("photonE", pho_e)
tree.SetBranchAddress("photonPx", pho_px)
tree.SetBranchAddress("photonPy", pho_py)
tree.SetBranchAddress("photonPz", pho_pz)
tree.SetBranchAddress("genPhotonE", genpho_e)
tree.SetBranchAddress("genPhotonPx", genpho_px)
tree.SetBranchAddress("genPhotonPy", genpho_py)
tree.SetBranchAddress("genPhotonPz", genpho_pz)
tree.SetBranchAddress("genPi0E", genpi0_e)
tree.SetBranchAddress("genPi0M", genpi0_m)

# Constants for pi0 mass selection
PI0_MASS = 0.135  # GeV
MASS_WINDOW = 0.05  # 50 MeV mass tolerance for π⁰

# Scan all gen photons and pi0 candidates to find min/max eta
photon_eta_min = 1e9
photon_eta_max = -1e9
pi0_eta_min = 1e9
pi0_eta_max = -1e9

pt_bins = array('d', [0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20])

for i_event in range(tree.GetEntries()):
    tree.GetEntry(i_event)
    gen_photons = [ROOT.TLorentzVector(genpho_px[j], genpho_py[j], genpho_pz[j], genpho_e[j]) for j in range(genpho_e.size())]
    for g in gen_photons:
        eta = g.Eta()
        if eta < photon_eta_min:
            photon_eta_min = eta
        if eta > photon_eta_max:
            photon_eta_max = eta
    # pi0 candidates
    for i in range(len(gen_photons)):
        for j in range(i+1, len(gen_photons)):
            p1 = gen_photons[i]
            p2 = gen_photons[j]
            pair = p1 + p2
            if abs(pair.M() - PI0_MASS) > MASS_WINDOW:
                continue
            eta = pair.Eta()
            if eta < pi0_eta_min:
                pi0_eta_min = eta
            if eta > pi0_eta_max:
                pi0_eta_max = eta

# Use dynamic binning for eta
def make_eta_bins(eta_min, eta_max, n_bins=20):
    # Use slightly padded range
    pad = 0.05 * (eta_max - eta_min)
    eta_min -= pad
    eta_max += pad
    return array('d', [eta_min + i*(eta_max-eta_min)/n_bins for i in range(n_bins+1)])

photon_eta_bins = make_eta_bins(photon_eta_min, photon_eta_max)
pi0_eta_bins = make_eta_bins(pi0_eta_min, pi0_eta_max)

# Create denominator and numerator histograms for efficiencies
photon_pt_den = ROOT.TH1D("photon_pt_den", "Photon pT;p_{T} [GeV];Efficiency", len(pt_bins)-1, pt_bins)
photon_pt_num = ROOT.TH1D("photon_pt_num", "Photon pT;p_{T} [GeV];Efficiency", len(pt_bins)-1, pt_bins)
photon_eta_den = ROOT.TH1D("photon_eta_den", "Photon #eta;#eta;Efficiency", len(photon_eta_bins)-1, photon_eta_bins)
photon_eta_num = ROOT.TH1D("photon_eta_num", "Photon #eta;#eta;Efficiency", len(photon_eta_bins)-1, photon_eta_bins)

pi0_pt_den = ROOT.TH1D("pi0_pt_den", "π^{0} pT;p_{T} [GeV];Efficiency", len(pt_bins)-1, pt_bins)
pi0_pt_num = ROOT.TH1D("pi0_pt_num", "π^{0} pT;p_{T} [GeV];Efficiency", len(pt_bins)-1, pt_bins)
pi0_eta_den = ROOT.TH1D("pi0_eta_den", "π^{0} #eta;#eta;Efficiency", len(pi0_eta_bins)-1, pi0_eta_bins)
pi0_eta_num = ROOT.TH1D("pi0_eta_num", "π^{0} #eta;#eta;Efficiency", len(pi0_eta_bins)-1, pi0_eta_bins)
# Constants:
PI0_MASS = 0.135  # GeV
MASS_WINDOW = 0.05  # 50 MeV mass tolerance for π⁰

reco_theta = []
n_matched_photons = 0
n_matched_pi0 = 0
# First, determine min and max theta of reconstructed photons
for i_event in range(tree.GetEntries()):
    tree.GetEntry(i_event)

    if len(genpho_e) == 0 or len(pho_e) == 0 or len(genpi0_m) == 0:
        continue

    # Construct TLorentzVectors
    reco_photons = [ROOT.TLorentzVector(pho_px[j], pho_py[j], pho_pz[j], pho_e[j]) for j in range(pho_e.size())]
    gen_photons = [ROOT.TLorentzVector(genpho_px[j], genpho_py[j], genpho_pz[j], genpho_e[j]) for j in range(genpho_e.size())]

    for reco in reco_photons:
        theta = reco.Theta()
        reco_theta.append(theta)

min_theta = min(reco_theta)
max_theta = max(reco_theta)
theta_cut_failed = 0
theta_cut_passed = 0
# Loop over events
for i_event in range(tree.GetEntries()):
    tree.GetEntry(i_event)

    if len(genpho_e) == 0 or len(genpi0_m) == 0:
        continue

    # Construct TLorentzVectors
    reco_photons = [ROOT.TLorentzVector(pho_px[j], pho_py[j], pho_pz[j], pho_e[j]) for j in range(pho_e.size())]
    gen_photons = [ROOT.TLorentzVector(genpho_px[j], genpho_py[j], genpho_pz[j], genpho_e[j]) for j in range(genpho_e.size())]        

    used_gen_indices = set()
    used_pi0_indices = set()

    # Pair gen photons with π⁰ candidates
    for i in range(len(gen_photons)):
        if i in used_gen_indices:
            continue
        # Energy cut for gen photon i
        if gen_photons[i].E() < 0.2:
            continue
        for j in range(i + 1, len(gen_photons)):
            # Energy cut for gen photon j
            if gen_photons[j].E() < 0.2:
                continue
            theta1 = gen_photons[i].Theta()
            theta2 = gen_photons[j].Theta()

            # Skip if either photon is outside detector acceptance
            if theta1 < min_theta or theta1 > max_theta:
                continue
            if theta2 < min_theta or theta2 > max_theta:
                continue

            if j in used_gen_indices:
                continue

            p1 = gen_photons[i]
            p2 = gen_photons[j]
            pair = p1 + p2  # Combined π0 candidate
            pair_mass = pair.M()

            if abs(pair_mass - PI0_MASS) > MASS_WINDOW:
                continue

            # Match pair to closest gen π⁰ by mass
            min_mass_diff = float('inf')
            best_pi0_idx = -1
            for k in range(len(genpi0_m)):
                if k in used_pi0_indices:
                    continue
                mass_diff = abs(pair_mass - genpi0_m[k])
                if mass_diff < min_mass_diff:
                    min_mass_diff = mass_diff
                    best_pi0_idx = k

            if best_pi0_idx >= 0:
                used_pi0_indices.add(best_pi0_idx)
                used_gen_indices.update([i, j])

                # Fill denominators for photons and pi0
                photon_pt_den.Fill(p1.Pt())
                photon_pt_den.Fill(p2.Pt())
                photon_eta_den.Fill(p1.Eta())
                photon_eta_den.Fill(p2.Eta())
                
                pi0_pt_den.Fill(pair.Pt())
                pi0_eta_den.Fill(pair.Eta())

                # Match each gen photon to reco photon
                used_reco_indices = set()
                matched_photons = []
                for gen_photon in [p1, p2]:
                    best_dr = float('inf')
                    best_reco_idx = -1
                    for idx, reco_photon in enumerate(reco_photons):
                        if idx in used_reco_indices:
                            continue
                        dr = reco_photon.DeltaR(gen_photon)
                        if dr < best_dr:
                            best_dr = dr
                            best_reco_idx = idx
                    if best_dr < 0.04 and best_reco_idx != -1:
                        matched_photons.append(gen_photon)
                        used_reco_indices.add(best_reco_idx)
                
                # Fill numerators for matched photons
                for gen_photon in matched_photons:
                    photon_pt_num.Fill(gen_photon.Pt())
                    photon_eta_num.Fill(gen_photon.Eta())
                
                # Fill π0 numerator only if both photons matched
                if len(matched_photons) == 2:
                    pi0_pt_num.Fill(pair.Pt())
                    pi0_eta_num.Fill(pair.Eta())
# Create TEfficiency objects for plotting
eff_photon_pt = ROOT.TEfficiency(photon_pt_num, photon_pt_den)
eff_photon_eta = ROOT.TEfficiency(photon_eta_num, photon_eta_den)
eff_pi0_pt = ROOT.TEfficiency(pi0_pt_num, pi0_pt_den)
eff_pi0_eta = ROOT.TEfficiency(pi0_eta_num, pi0_eta_den)


ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptStat(0)

canvas_pt = ROOT.TCanvas("canvas_pt", "Efficiency vs p_{T}", 1200, 900)
canvas_eta = ROOT.TCanvas("canvas_eta", "Efficiency vs #eta", 1200, 900)


eff_photon_pt.SetLineColor(ROOT.kBlue)
eff_photon_pt.SetMarkerColor(ROOT.kBlue)
eff_photon_pt.SetMarkerStyle(20)
eff_photon_pt.SetMarkerSize(1.2)

eff_pi0_pt.SetLineColor(ROOT.kRed)
eff_pi0_pt.SetMarkerColor(ROOT.kRed)
eff_pi0_pt.SetMarkerStyle(21)
eff_pi0_pt.SetMarkerSize(1.2)

eff_photon_eta.SetLineColor(ROOT.kBlue)
eff_photon_eta.SetMarkerColor(ROOT.kBlue)
eff_photon_eta.SetMarkerStyle(20)
eff_photon_eta.SetMarkerSize(1.2)

eff_pi0_eta.SetLineColor(ROOT.kRed)
eff_pi0_eta.SetMarkerColor(ROOT.kRed)
eff_pi0_eta.SetMarkerStyle(21)
eff_pi0_eta.SetMarkerSize(1.2)

# Find min and max efficiency values for pT plots
photon_pt_graph = eff_photon_pt.CreateGraph()
pi0_pt_graph = eff_pi0_pt.CreateGraph()

min_eff = float('inf')
max_eff = float('-inf')

# Check photon efficiencies
for i in range(photon_pt_graph.GetN()):
    eff = photon_pt_graph.GetY()[i]
    if eff > 0:  # Ignore zero efficiencies
        min_eff = min(min_eff, eff)
        max_eff = max(max_eff, eff)

# Check pi0 efficiencies
for i in range(pi0_pt_graph.GetN()):
    eff = pi0_pt_graph.GetY()[i]
    if eff > 0:  # Ignore zero efficiencies
        min_eff = min(min_eff, eff)
        max_eff = max(max_eff, eff)

# Add 10% padding to the range
padding = (max_eff - min_eff) * 0.1
y_min = max(0, min_eff - padding)
y_max = min(1, max_eff + padding)

# Draw pT efficiencies
canvas_pt.cd()
canvas_pt.SetLogy(0)  # Make sure we're in linear scale
eff_photon_pt.SetTitle("Reconstruction Efficiency vs p_{T}")

# Create a histframe to control the range
hframe = ROOT.TH1F("hframe", "", 100, pt_bins[0], pt_bins[-1])
hframe.GetYaxis().SetRangeUser(y_min, y_max)  # Set the dynamic y-axis range
hframe.SetTitle("Reconstruction Efficiency vs p_{T}")
hframe.GetXaxis().SetTitle("p_{T} [GeV]")
hframe.GetYaxis().SetTitle("Efficiency")
hframe.Draw()

# Draw efficiencies on top
eff_photon_pt.Draw("SAME")
eff_pi0_pt.Draw("SAME")

# Create legend for pT plot
legend_pt = ROOT.TLegend(0.65, 0.15, 0.85, 0.35)
legend_pt.SetBorderSize(0)
legend_pt.SetFillStyle(0)
legend_pt.AddEntry(eff_photon_pt, "Single Photon", "lep")
legend_pt.AddEntry(eff_pi0_pt, "#pi^{0}", "lep")
legend_pt.Draw()

# Draw eta efficiencies
canvas_eta.cd()
eff_photon_eta.SetTitle("Reconstruction Efficiency vs #eta")
# Set y-axis range for eta plot
eta_graph = eff_photon_eta.CreateGraph()
pi0_eta_graph = eff_pi0_eta.CreateGraph()

min_eff = float('inf')
max_eff = float('-inf')

# Check photon efficiencies
for i in range(eta_graph.GetN()):
    eff = eta_graph.GetY()[i]
    if eff > 0:  # Ignore zero efficiencies
        min_eff = min(min_eff, eff)
        max_eff = max(max_eff, eff)

# Check pi0 efficiencies
for i in range(pi0_eta_graph.GetN()):
    eff = pi0_eta_graph.GetY()[i]
    if eff > 0:  # Ignore zero efficiencies
        min_eff = min(min_eff, eff)
        max_eff = max(max_eff, eff)

# Add 10% padding to the range
padding = (max_eff - min_eff) * 0.1
y_min = max(0, min_eff - padding)
y_max = min(1, max_eff + padding)

hframe_eta = ROOT.TH1F("hframe_eta", "", 100, photon_eta_bins[0], photon_eta_bins[-1])
hframe_eta.GetYaxis().SetRangeUser(y_min, y_max)  # Set dynamic y-axis range
hframe_eta.SetTitle("Reconstruction Efficiency vs #eta")
hframe_eta.GetXaxis().SetTitle("#eta")
hframe_eta.GetYaxis().SetTitle("Efficiency")
hframe_eta.Draw()
eff_photon_eta.Draw()
eff_pi0_eta.Draw("SAME")

# Create legend for eta plot
legend_eta = ROOT.TLegend(0.65, 0.15, 0.85, 0.35)
legend_eta.SetBorderSize(0)
legend_eta.SetFillStyle(0)
legend_eta.AddEntry(eff_photon_eta, "Single Photon", "lep")
legend_eta.AddEntry(eff_pi0_eta, "#pi^{0}", "lep")
legend_eta.Draw()

canvas_pt.Update()
canvas_eta.Update()
canvas_pt.SaveAs("efficiency_pt.png")
canvas_eta.SaveAs("efficiency_eta.png")

print("Plots saved as efficiency_pt.png and efficiency_eta.png")

