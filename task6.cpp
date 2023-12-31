#include <iostream>
#include <TMath.h>
#include <TF1.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include <TVector3.h>
#include <TH1D.h>


const Double_t full_energy_MeV = 1020;
const Double_t e_mass_MeV = 0.51099895;
const Double_t K_mass_MeV = 497.611;
const Double_t pi_mass_MeV = 139.57039;
const Int_t N = 1000;
const Double_t pi_ctau = 7.8045; // метры
const Double_t K_ctau = 0.026; // метры
const Double_t detector_R = 0.3; // метры
const Double_t detector_L = 0.5; // метры

Double_t Get_px_from_E(Double_t E, Double_t m) 
{
    return E * TMath :: Power(1 - TMath :: Power(m / E, 2),  0.5);
}

Double_t Get_gamma_from_E(Double_t E, Double_t m)
{
    return E / (m);
}



void task6()
{
    auto c = new TCanvas("task6_pi", "task6_pi", 1000, 1000);
    auto c1 = new TCanvas("task6_Ks", "task6_Ks", 1000, 1000);
    c->Divide(4, 1);
    c1->Divide(2, 1);
    auto Ks_polar = new TH1D("Ks_polar", "Ks_polar", 50, 0, TMath :: Pi());
    auto Ks_azimuthal = new TH1D("Ks_azimuthal", "Ks_azimuthal", 50, -TMath::Pi(), TMath :: Pi());



    auto pi_polar_Ks_hist = new TH1D("pi_polar_Ks", "pi_polar_Ks", 50, 0, TMath :: Pi());
    auto Ks_flyL = new TH1D("Ks_uflyL", "Ks_flyL;m;Entries", 50, 0, 0);
    auto pi_polar_ls_hist = new TH1D("pi_polar_ls", "pi_polar_ls", 50, 0, TMath :: Pi());
    auto pi_azimuthal_ls_hist = new TH1D("pi_azimuthal_ls", "pi_azimuthal_ls", 50, -TMath::Pi(), TMath :: Pi());
    Int_t cnt = 0;


    TF1 *func = new TF1("f", "(sin(x))^3", 0, TMath::Pi());
    TF1 *ctau_func = new TF1("ctau_f", "exp(-x/[0])", 0, 50);



    for (Int_t i = 0; i < N; i ++)
    {
        std::cout << i << "'s iteration has started" << std::endl;
        
        // генерация Ks (Kl) в ЛСО 

        Double_t theta_K = func->GetRandom(); 
        Double_t phi_K = gRandom->Rndm() * 2 * TMath :: Pi(); 

        Double_t Ks_px = Get_px_from_E(full_energy_MeV / 2., K_mass_MeV); 
        TLorentzVector  Ks(Ks_px, 0, 0, full_energy_MeV / 2.);
        ctau_func->SetParameter(0, K_ctau);
        TLorentzVector  Ks_decay(0, 0, 0, ctau_func->GetRandom());
       
        


        // поворот 4-х вектора Ks на сгенерированное направление

        Ks.SetPhi(phi_K);
        Ks.SetTheta(theta_K);
        
        Ks_polar->Fill(Ks.Theta());
        Ks_azimuthal->Fill(Ks.Phi());

        std :: cout << "Ks in lab frame before boosts: {" << Ks[0] << ", " << Ks[1] << ", " << Ks[2] << ", " << Ks[3] << "}" << std :: endl;

        TLorentzVector  Kl = Ks; 
        Kl.RotateZ(TMath::Pi());
        

        // переход в систему отчета Ks
        Ks_decay.Boost(Ks.BoostVector());
        TLorentzRotation from_lab_to_Ks;
        from_lab_to_Ks.Boost(-Ks.BoostVector());
        Ks = from_lab_to_Ks * Ks;
        Double_t new_full_energy_MeV = Ks.E();

        std :: cout << "Ks in Ks frame: {" << Ks[0] << ", " << Ks[1] << ", " << Ks[2] << ", " << Ks[3] << "}" << std :: endl;

        // генерация pi+ (-) в системе Ks
        Double_t theta_pi = TMath :: ACos(gRandom->Rndm() * 2 - 1);
        Double_t phi_pi = gRandom->Rndm() * 2 * TMath :: Pi(); 

        Double_t pi_plus_px = Get_px_from_E(new_full_energy_MeV / 2., pi_mass_MeV); 
        TLorentzVector  pi_plus(pi_plus_px, 0, 0, new_full_energy_MeV / 2.);

        
        

        // время жизни
        ctau_func->SetParameter(0, pi_ctau);
        TLorentzVector  pi_decay_plus(0, 0, 0, ctau_func->GetRandom()); 
        TLorentzVector  pi_decay_minus(0, 0, 0, ctau_func->GetRandom());
        
        pi_decay_plus.Boost(pi_plus.BoostVector());
        pi_decay_minus.Boost(pi_plus.BoostVector());
    


        pi_decay_plus.SetPhi(phi_pi);
        pi_decay_plus.SetTheta(theta_pi);
        pi_decay_minus.SetPhi(phi_pi);
        pi_decay_minus.SetTheta(-theta_pi);
        pi_plus.SetPhi(phi_pi);
        pi_plus.SetTheta(theta_pi);

        //Заполнение гистограммы в системе Ks 
        pi_polar_Ks_hist->Fill(pi_plus.Theta());

        
       TLorentzVector pi_minus;
       pi_minus[0] = - pi_plus[0];
       pi_minus[1] = - pi_plus[1];
       pi_minus[2] = - pi_plus[2];
       pi_minus[3] =  pi_plus[3];
       std :: cout << "Pi_minus Pi_plus summ in Ks frame: {" << pi_minus[0] + pi_plus[0] << ", " << pi_minus[1] + pi_plus[1] << ", " << pi_minus[2] + pi_plus[2] << ", " << pi_minus[3] + pi_plus[3] << "}"<< std::endl;



        // перевод векторов в ЛСО
        TLorentzRotation from_Ks_to_lab = from_lab_to_Ks.Inverse();
        pi_minus = from_Ks_to_lab * pi_minus;
        pi_plus = from_Ks_to_lab * pi_plus;
        pi_decay_plus = from_Ks_to_lab * pi_decay_plus;
        pi_decay_minus = from_Ks_to_lab * pi_decay_minus;
        Ks = from_Ks_to_lab * Ks;


        //Заполнение гистограммы в системе ls
        pi_polar_ls_hist->Fill(pi_plus.Theta());
        pi_azimuthal_ls_hist->Fill(pi_plus.Phi());


        Ks_flyL->Fill(TMath :: Power(Ks_decay.X() * Ks_decay.X() + Ks_decay.Y() * Ks_decay.Y() + Ks_decay.Z() * Ks_decay.Z(), 0.5));

        std :: cout << "Ks in lab frame after boosts: {" << Ks[0] << ", " << Ks[1] << ", " << Ks[2] << ", " << Ks[3] << "}" << std :: endl;
        std :: cout << "Pi_minus Pi_plus summ in lab frame: {" << pi_minus[0] + pi_plus[0] << ", " << pi_minus[1] + pi_plus[1] << ", " << pi_minus[2] + pi_plus[2] << ", " << pi_minus[3] + pi_plus[3] << "}"<< std::endl;

        Double_t pi_plus_P_transverse = TMath::Sqrt(pi_plus.Px() * pi_plus.Px() + pi_plus.Py() * pi_plus.Py());
        Double_t pi_minus_P_transverse = TMath::Sqrt(pi_minus.Px() * pi_minus.Px() + pi_minus.Py() * pi_minus.Py());
        Double_t Ks_decay_transverse =  TMath::Sqrt(Ks_decay.X() * Ks_decay.X() + Ks_decay.Y() * Ks_decay.Y());
        Double_t pi_decay_plus_transverse = TMath::Sqrt(pi_decay_plus.X() * pi_decay_plus.X() + pi_decay_plus.Y() * pi_decay_plus.Y());
        Double_t pi_decay_minus_transverse = TMath::Sqrt(pi_decay_minus.X() * pi_decay_minus.X() + pi_decay_minus.Y() * pi_decay_minus.Y());

        if (((pi_plus_P_transverse > 40) && pi_plus.Theta() > TMath::ATan(detector_R / (detector_L - Ks_decay_transverse - pi_decay_plus_transverse))) &&
        ((pi_minus_P_transverse > 40) && pi_minus.Theta() > TMath::ATan(detector_R / (detector_L - Ks_decay_transverse - pi_decay_minus_transverse)))) cnt ++; 
        

    }
    c->cd(1);
    pi_polar_Ks_hist->Draw();
    c->cd(2);
    Ks_flyL->Draw();
    c->cd(3);
    pi_polar_ls_hist->Draw();
    c->cd(4);
    pi_azimuthal_ls_hist->Draw();

    c1->cd(1);
    Ks_azimuthal->Draw();
    c1->cd(2);
    Ks_polar->Draw();

    Double_t bin_p = Double_t(cnt) / N;
    std :: cout << "efficency: " << Double_t(cnt) / N <<  " +- " << TMath::Sqrt(bin_p * (1 - bin_p)) << std :: endl;
    

}
