#include "Mt2Calculator.h"

#include <cmath>
#include <algorithm>
#include <vector>

#include "TLorentzVector.h"

namespace Mt2Calculator{
	Calculator::Calculator() : metResM(TMatrixD(2,2))
	{
		neutrino_test_mass=0.0L;
		W_test_mass=80.0L;//COME BACK AND CHANGE THIS BACK LATER
		jet1=TLorentzVector(0.0, 0.0, 0.0, 0.0);
		jet2=TLorentzVector(0.0, 0.0, 0.0, 0.0);
		lep1=TLorentzVector(0.0, 0.0, 0.0, 0.0);
		lep2=TLorentzVector(0.0, 0.0, 0.0, 0.0);
		met=TLorentzVector(0.0, 0.0, 0.0, 0.0);
		jet1Uncor=TLorentzVector(0.0, 0.0, 0.0, 0.0);
		jet2Uncor=TLorentzVector(0.0, 0.0, 0.0, 0.0);
		jet1Res[0]=0.0L;
		jet1Res[1]=0.0L;
		jet1Res[2]=0.0L;
		jet1UncorRes[0]=0.0L;
		jet1UncorRes[1]=0.0L;
		jet1UncorRes[2]=0.0L;
		jet2Res[0]=0.0L;
		jet2Res[1]=0.0L;
		jet2Res[2]=0.0;
		jet2UncorRes[0]=0.0L;
		jet2UncorRes[1]=0.0L;
		jet2UncorRes[2]=0.0L;
		lep1Res[0]=0.0L;
		lep1Res[1]=0.0L;
		lep1Res[2]=0.0L;
		lep2Res[0]=0.0L;
		lep2Res[1]=0.0L;
		lep2Res[2]=0.0L;
		metRes[0][0]=0.0L;
		metRes[0][1]=0.0L;
		metRes[1][0]=0.0L;
		metRes[1][1]=0.0L;
	}

	void Calculator::MatrixMultiply(const long double first[][2],
			const long double second[][2], long double product[][2]){
		product[0][0]=first[0][0]*second[0][0]+first[0][1]*second[1][0];
		product[0][1]=first[0][0]*second[0][1]+first[0][1]*second[1][1];
		product[1][0]=first[1][0]*second[0][0]+first[1][1]*second[1][0];
		product[1][1]=first[1][0]*second[0][1]+first[1][1]*second[1][1];
	}

	void Calculator::SetParticles(const TLorentzVector &jet1in, const TLorentzVector &jet2in,
			const TLorentzVector &lep1in, const TLorentzVector &lep2in,
			const TLorentzVector &metin){
		jet1=jet1in;
		jet2=jet2in;
		lep1=lep1in;
		lep2=lep2in;
		met=metin;
	}

	void Calculator::SetUncorrectedJets(const TLorentzVector &jet1in,
			const TLorentzVector &jet2in){
		jet1Uncor=jet1in;
		jet2Uncor=jet2in;
	}

	void Calculator::SetNeutrinoTestMass(const double mcin){
		neutrino_test_mass=(long double)(mcin);
	}

	void Calculator::SetWTestMass(const double mcin){
		W_test_mass=(long double)(mcin);
	}

	void Calculator::SetResolutions(const double jet1ResIn[], const double jet1ResUncorIn[],
			const double jet2ResIn[], const double jet2ResUncorIn[],
			const double lep1ResIn[], const double lep2ResIn[],
			const TMatrixD &metResIn){
		jet1Res[0]=jet1ResIn[0];
		jet1Res[1]=jet1ResIn[1];
		jet1Res[2]=jet1ResIn[2];
		jet1UncorRes[0]=jet1ResUncorIn[0];
		jet1UncorRes[1]=jet1ResUncorIn[1];
		jet1UncorRes[2]=jet1ResUncorIn[2];
		jet2Res[0]=jet2ResIn[0];
		jet2Res[1]=jet2ResIn[1];
		jet2Res[2]=jet2ResIn[2];
		jet2UncorRes[0]=jet2ResUncorIn[0];
		jet2UncorRes[1]=jet2ResUncorIn[1];
		jet2UncorRes[2]=jet2ResUncorIn[2];
		lep1Res[0]=lep1ResIn[0];
		lep1Res[1]=lep1ResIn[1];
		lep1Res[2]=lep1ResIn[2];
		lep2Res[0]=lep2ResIn[0];
		lep2Res[1]=lep2ResIn[1];
		lep2Res[2]=lep2ResIn[2];
		metRes[0][0]=(long double)(metResIn(0,0));
		metRes[0][1]=(long double)(metResIn(0,1));
		metRes[1][0]=(long double)(metResIn(1,0));
		metRes[1][1]=(long double)(metResIn(1,1));
		metResM = metResIn;
	}


	double Calculator::GetUpstream210PhiResolution() {
		// We're not going to try to remove the lepton resolutions
		// Totate the upstream covariance matrix so that the the pt is along the x-axis
		TMatrixD upstreamCovMatrix(2,2);
		upstreamCovMatrix = metResM;
		TLorentzVector p4Up = -lep1-lep2-met;
		TMatrixD R(2,2);
		R(0,0) = cos(-p4Up.Phi());
		R(0,1) = -sin(-p4Up.Phi());
		R(1,0) = sin(-p4Up.Phi());
		R(1,1) = cos(-p4Up.Phi());
		TMatrixD RInv = R;
		RInv.Invert();
		upstreamCovMatrix = R*upstreamCovMatrix*RInv;
		return sqrt(upstreamCovMatrix(1,1))/p4Up.Pt();
	}

	double Calculator::GetUpstreamPhiResolution(){
		// covariance matrix for the two jets
		TMatrixD V1 = GetObjCovarianceMatrix(jet1Uncor, jet1UncorRes[0], jet1UncorRes[1]*jet1Uncor.Pt());
		TMatrixD V2 = GetObjCovarianceMatrix(jet2Uncor, jet2UncorRes[0], jet2UncorRes[1]*jet2Uncor.Pt());
		// remove contributions of two jets from the covariance matrix
		TMatrixD upstreamCovMatrix(2,2);
		upstreamCovMatrix = metResM-V1-V2;
		// For now, ignore contribution of leptons to this matrix
		
		// Now rotate the upstream covariance matrix so that the the pt is along the x-axis
		TLorentzVector p4Up = -jet1Uncor-jet2Uncor-lep1-lep2-met;
		TMatrixD R(2,2);
		R(0,0) = cos(-p4Up.Phi());
		R(0,1) = -sin(-p4Up.Phi());
		R(1,0) = sin(-p4Up.Phi());
		R(1,1) = cos(-p4Up.Phi());
		TMatrixD RInv = R;
		RInv.Invert();
		upstreamCovMatrix = R*upstreamCovMatrix*RInv;
		return sqrt(upstreamCovMatrix(1,1))/p4Up.Pt();
	}

	TMatrixD Calculator::GetObjCovarianceMatrix(TLorentzVector &p, double pt_res, double phi_res) {
		TMatrixD m(2,2);
		double cosphi = cos(p.Phi());
		double sinphi = sin(p.Phi());
		double Ctt = pt_res*pt_res;
		double Cff = phi_res*phi_res;

		m(0,0) = Ctt*cosphi*cosphi + Cff*sinphi*sinphi;
		m(0,1) = cosphi*sinphi*(Ctt-Cff);
		m(1,0) = m(0,1);
		m(1,1) = Cff*cosphi*cosphi + Ctt*sinphi*sinphi;

		return m;
	}

	TLorentzVector Calculator::W1a(){
		return jet1+lep1;
	}

	TLorentzVector Calculator::W2a(){
		return jet2+lep2;
	}

	TLorentzVector Calculator::W1b(){
		return jet2+lep1;
	}

	TLorentzVector Calculator::W2b(){
		return jet1+lep2;
	}

	double Calculator::GetMt2(const int p, const int c){
		if (p==1 && c==0){
			return CalcMt2Numerically(lep1, lep2, met, neutrino_test_mass);
		}else if (p==2 && c==0){
			const double a=CalcMt2Numerically(W1a(), W2a(), met, neutrino_test_mass);
			const double b=CalcMt2Numerically(W1b(), W2b(), met, neutrino_test_mass);
			return std::min(a,b);
		}else if (p==2 && c==1){
			return CalcMt2Numerically(jet1, jet2, lep1+lep2+met, W_test_mass);
		}else{
			return -1.0;
		}
	}

	double Calculator::GetMt2Perp(const int p, const int c){
		if (p==1 && c==0){
			return CalcMt2Perp(lep1, lep2, met, neutrino_test_mass);
		}else if (p==2 && c==0){
			const double a=CalcMt2Perp(W1a(), W2a(), met, neutrino_test_mass);
			const double b=CalcMt2Perp(W1b(), W2b(), met, neutrino_test_mass);
			return std::min(a,b);
		}else if (p==2 && c==1){
			return CalcMt2Perp(jet1, jet2, lep1+lep2+met, W_test_mass);
		}else{
			return -1.0;
		}
	}

	double Calculator::GetMctPerp(const int p, const int c, const double phi){
		if (p==1 && c==0){
			return CalcMctPerp(lep1, lep2, (long double)(phi));
		}else if (p==2 && c==0){
			const double a=CalcMctPerp(W1a(), W2a(), (long double)(phi));
			const double b=CalcMctPerp(W1b(), W2b(), (long double)(phi));
			return std::min(a,b);
		}else if (p==2 && c==1){
			return CalcMctPerp(jet1, jet2, (long double)(phi));
		}else{
			return -1.0;
		}
	}

	double Calculator::GetMct(const int p, const int c){
		if (p==1 && c==0){
			return CalcMct(lep1, lep2);
		}else if (p==2 && c==0){
			const double a=CalcMct(W1a(), W2a());
			const double b=CalcMct(W1b(), W2b());
			return std::min(a,b);
		}else if (p==2 && c==1){
			return CalcMct(jet1, jet2);
		}else{
			return -1.0;
		}
	}

	double Calculator::GetMctPerp(const int p, const int c){
		if (p==1 && c==0){
			return CalcMctPerp(-(lep1+lep2+met), lep1, lep2);
		}else if (p==2 && c==0){
			TLorentzVector upstream=-(jet1+jet2+lep1+lep2+met);
			const double a=CalcMctPerp(upstream, W1a(), W2a());
			const double b=CalcMctPerp(upstream, W1b(), W2b());
			return std::min(a,b);
		}else if (p==2 && c==1){
			return CalcMctPerp(-(jet1+jet2+lep1+lep2+met), jet1, jet2);
		}else{
			return -1.0;
		}
	}

	std::vector<double> Calculator::GetAllMt2_220(){
		const double mt2_220a=CalcMt2Numerically(W1a(), W2a(), met, neutrino_test_mass);
		const double mt2_220b=CalcMt2Numerically(W1b(), W2b(), met, neutrino_test_mass);
		std::vector<double> output;
		output.clear();
		output.push_back(mt2_220a);
		output.push_back(mt2_220b);
		return output;		
	}

	std::vector<double> Calculator::GetAllMt2Perp_220(){
		const double mt2_220a=CalcMt2Perp(W1a(), W2a(), met, neutrino_test_mass);
		const double mt2_220b=CalcMt2Perp(W1b(), W2b(), met, neutrino_test_mass);
		std::vector<double> output;
		output.clear();
		output.push_back(mt2_220a);
		output.push_back(mt2_220b);
		return output;
	}

	std::vector<double> Calculator::GetAllMct_220(){
		const double mct_220a=CalcMct(W1a(), W2a());
		const double mct_220b=CalcMct(W1b(), W2b());
		std::vector<double> output;
		output.clear();
		output.push_back(mct_220a);
		output.push_back(mct_220b);
		return output;		
	}

	std::vector<double> Calculator::GetAllMctPerp_220(){
		const TLorentzVector upstream=-(jet1+jet2+lep1+lep2+met);
		const double mct_220a=CalcMctPerp(upstream, W1a(), W2a());
		const double mct_220b=CalcMctPerp(upstream, W1b(), W2b());
		std::vector<double> output;
		output.clear();
		output.push_back(mct_220a);
		output.push_back(mct_220b);
		return output;
	}

	std::vector<double> Calculator::GetBlInvariantMasses(){
		const long double A=std::max(W1a().M(), W2a().M());
		const long double a=std::min(W1a().M(), W2a().M());
		const long double B=std::max(W1b().M(), W2b().M());
		const long double b=std::min(W1b().M(), W2b().M());

		std::vector<double> output;
		output.clear();

		if (A>=a && a>=B && B>=b){
			output.push_back(B);
			output.push_back(b);
		}else if (A>=B && B>=a && a>=b){
			output.push_back(B);
			output.push_back(a);
			output.push_back(b);
		}else if (A>=B && B>=b && b>=a){
			output.push_back(B);
			output.push_back(a);
			output.push_back(b);
		}else if (B>=b && b>=A && A>=a){
			output.push_back(A);
			output.push_back(a);
		}else if (B>=A && A>=b && b>=a){
			output.push_back(A);
			output.push_back(a);
			output.push_back(b);
		}else{//B>=A && A>=a && a>=b
			output.push_back(A);
			output.push_back(a);
			output.push_back(b);
		}

		return output;
	}

	std::vector<double> Calculator::GetBMasses(){
		const long double A=std::max(W1a().M(), W2a().M());
		const long double a=std::min(W1a().M(), W2a().M());
		const long double B=std::max(W1b().M(), W2b().M());
		const long double b=std::min(W1b().M(), W2b().M());
		
		std::vector<double> output;
		output.clear();

		long double Amass=0, amass=0, Bmass=0, bmass=0;
		if(W1a().M() > W2a().M()){
			Amass=jet1.M();
			amass=jet2.M();
		}
		if(W1a().M() <= W2a().M()){
			Amass=jet2.M();
			amass=jet1.M();
		}
		if(W1b().M() > W2b().M()){
			Bmass=jet2.M();
			bmass=jet1.M();
		}
		if(W1b().M() <= W2b().M()){
			Bmass=jet1.M();
			bmass=jet2.M();
		}

		if (A>=a && a>=B && B>=b){
			output.push_back(Bmass);
			output.push_back(bmass);
		}else if (A>=B && B>=a && a>=b){
			output.push_back(Bmass);
			output.push_back(amass);
			output.push_back(bmass);
		}else if (A>=B && B>=b && b>=a){
			output.push_back(Bmass);
			output.push_back(amass);
			output.push_back(bmass);
		}else if (B>=b && b>=A && A>=a){
			output.push_back(Amass);
			output.push_back(amass);
		}else if (B>=A && A>=b && b>=a){
			output.push_back(Amass);
			output.push_back(amass);
			output.push_back(bmass);
		}else{//B>=A && A>=a && a>=b
			output.push_back(Amass);
			output.push_back(amass);
			output.push_back(bmass);
		}
		return output;
	}

	double Calculator::GetMctPerp221Resolution(){ 

		const TLorentzVector b1=jet1;
		const TLorentzVector b2=jet2;
		const long double m1 = b1.M();
		const long double m2 = b2.M();
		const long double pT1 = b1.Pt();
		const long double pT2 = b2.Pt();
		const TLorentzVector Up=jet1+jet2+lep1+lep2+met;
		const long double delta1 = b1.DeltaPhi(Up);
		const long double delta2 = b2.DeltaPhi(Up);
		const long double resPhiup=GetUpstreamPhiResolution();
		const long double resPT1 = jet1Res[0];
		const long double resPhi1 = jet1Res[1];
		const long double resPT2 = jet2Res[0];
		const long double resPhi2 = jet2Res[1];

		const long double ET1 = sqrt(m1*m1+pT1*pT1*sin(delta1)*sin(delta1));
		const long double ET2 = sqrt(m2*m2+pT2*pT2*sin(delta2)*sin(delta2));
		long double mct = m1*m1+m2*m2+2*(ET1*ET2+pT1*pT2*sin(delta1)*sin(delta2));
		mct = sqrt(mct);
		const long double resDel1 = sqrt(resPhi1*resPhi1+resPhiup*resPhiup);
		const long double resDel2 = sqrt(resPhi2*resPhi2+resPhiup*resPhiup);
		const long double dET1dpT1 = pT1*sin(delta1)*sin(delta1)/ET1;
		const long double dET2dpT2 = pT2*sin(delta2)*sin(delta2)/ET2;
		const long double dET1ddel1 = pT1*pT1*sin(delta1)*cos(delta1)/ET1;
		const long double dET2ddel2 = pT2*pT2*sin(delta2)*cos(delta2)/ET2;
		const long double dMCTdpt1 = (dET1dpT1*ET2+pT2*sin(delta1)*sin(delta2))/mct;
		const long double dMCTdpt2 = (dET2dpT2*ET1+pT1*sin(delta2)*sin(delta1))/mct;
		const long double dMCTddel1 = (dET1ddel1*ET2+pT1*pT2*sin(delta2)*cos(delta1))/mct;
		const long double dMCTddel2 = (dET2ddel2*ET1+pT2*pT1*sin(delta1)*cos(delta2))/mct;

		const long double out = dMCTdpt1*dMCTdpt1*resPT1*resPT1+dMCTdpt2*dMCTdpt2*resPT2*resPT2+dMCTddel1*dMCTddel1*resDel1*resDel1+dMCTddel2*dMCTddel2*resDel2*resDel2;

		return sqrt(out);
	}

	double Calculator::GetMt2Perp221Resolution(){
		const TLorentzVector &b1=jet1;
		const TLorentzVector &b2=jet2;
		const TLorentzVector &Up=jet1+jet2+lep1+lep2+met;
		const long double sigUp2=GetUpstreamPhiResolution();
		const long double MT2=GetMt2Perp(2,1);
		const long double sigPT = jet1Res[0];
		const long double sigTh = jet1Res[1];
		const long double sigPT2 = jet2Res[0];
		const long double sigTh2 = jet2Res[1];
		const long double mc=W_test_mass;

		long double alpha = 1.0L;	
		TVector3 hold1 = b1.Vect().Cross(Up.Vect());
		TVector3 hold2 = b2.Vect().Cross(Up.Vect());
		if (hold1.Dot(hold2) < 0.0L){
			alpha = -1.0L;
		}
		const long double pt1 = b1.Pt()*fabs(sin(b1.DeltaPhi(Up)));
		const long double pt2 = b2.Pt()*fabs(sin(b2.DeltaPhi(Up)));
		const long double ET1 = sqrt(b1.M2()+pt1*pt1);
		const long double ET2 = sqrt(b2.M2()+pt2*pt2);
		const long double AT = ET1*ET2+alpha*pt1*pt2;
		const long double msq = b1.M2()+b2.M2();
		const long double M = b1.M()*b2.M();
		const long double numerator = 4.0L*AT*AT*AT+AT*AT*(4.0L*msq+4.0L*mc*mc)+AT*(msq*msq+4.0L*mc*mc*msq)+4.0L*mc*mc*M*M;
		const long double denominator = (MT2*MT2-mc*mc-AT)*(2.0L*AT+msq)*(2.0L*AT+msq);
		const long double part1 = 0.5L*(numerator/denominator + 1.0L)/MT2;
		const long double parta = pt1/(ET1*b1.Pt())*(ET2*pt1+alpha*ET1*pt2)*sigPT;
		const long double partb = pt2/(ET2*b2.Pt())*(ET1*pt2+alpha*ET2*pt1)*sigPT2;
		const long double partc = pt1/(ET1*tan(b1.DeltaPhi(Up)))*(ET2*pt1+alpha*ET1*pt2)*sqrt(sigTh*sigTh+sigUp2*sigUp2);
		const long double partd = pt2/(ET2*tan(b2.DeltaPhi(Up)))*(ET1*pt2+alpha*ET2*pt1)*sqrt(sigTh2*sigTh2+sigUp2*sigUp2);
		return sqrt(part1*part1*(parta*parta+partb*partb+partc*partc+partd*partd));
	}

	double Calculator::ResIM(const TLorentzVector &b1, const TLorentzVector &l1, const long double bres[]){
		const long double sigPT = bres[0];
		const long double sigTh = bres[1];
		const long double sigEta = bres[2];

		TLorentzVector bl = b1+l1;
		TVector3 Ptl = TVector3(l1.Px(), l1.Py(),0.0);
		TVector3 Ptb = TVector3(b1.Px(), b1.Py(),0.0);
		TVector3 Pl = l1.Vect();	
		TVector3 Pb = b1.Vect();
		const long double betab = b1.Beta();
		const long double term1 = sigPT/Ptb.Mag()*(Pb.Dot(Pl)-Pl.Mag()*Pb.Mag()*betab);
		const long double term2 = Ptb.Mag()*Ptl.Mag()*sin(Ptb.Phi()-Ptl.Phi())*sigTh;
		const long double term3 = sigEta*betab*(l1.Pz()*b1.E()-Pl.Mag()*b1.Pz());

		const long double out = (term1*term1+term2*term2+term3*term3)/(bl.M2());
		return sqrt(out);
	}

	std::vector<double> Calculator::GetBlInvariantMassResolutions(){
		long double A, a, B, b, resA[3], resa[3], resB[3], resb[3];
		TLorentzVector bA, ba, bB, bb, lA, la, lB, lb;
		if (W1a().M() > W2a().M()){
			A=W1a().M();
			bA=jet1;
			lA=lep1;
			resA[0]=jet1Res[0];
			resA[1]=jet1Res[1];
			resA[2]=jet1Res[2];
			a=W2a().M();
			ba=jet2;
			la=lep2;
			resa[0]=jet2Res[0];
			resa[1]=jet2Res[1];
			resa[2]=jet2Res[2];
		}else{
			A=W2a().M();
			bA=jet2;
			lA=lep2;
			resA[0]=jet2Res[0];
			resA[1]=jet2Res[1];
			resA[2]=jet2Res[2];
			a=W1a().M();
			ba=jet1;
			la=lep1;
			resa[0]=jet1Res[0];
			resa[1]=jet1Res[1];
			resa[2]=jet1Res[2];
		}
		if (W1b().M() > W2b().M()){
			B=W1b().M();
			bB=jet2;
			lB=lep1;
			resB[0]=jet2Res[0];
			resB[1]=jet2Res[1];
			resB[2]=jet2Res[2];
			b=W2b().M();
			bb=jet1;
			lb=lep2;
			resb[0]=jet1Res[0];
			resb[1]=jet1Res[1];
			resb[2]=jet1Res[2];
		}else{
			B=W2b().M();
			bB=jet1;
			lB=lep2;
			resB[0]=jet1Res[0];
			resB[1]=jet1Res[1];
			resB[2]=jet1Res[2];
			b=W1b().M();
			bb=jet2;
			lb=lep1;
			resb[0]=jet2Res[0];
			resb[1]=jet2Res[1];
			resb[2]=jet2Res[2];
		}

		std::vector<double> output;
		output.clear();

		if (A>=a && a>=B && B>=b){
			output.push_back(ResIM(bB, lB, resB));
			output.push_back(ResIM(bb, lb, resb));
		}else if (A>=B && B>=a && a>=b){
			output.push_back(ResIM(bB, lB, resB));
			output.push_back(ResIM(ba, la, resa));
			output.push_back(ResIM(bb, lb, resb));
		}else if (A>=B && B>=b && b>=a){
			output.push_back(ResIM(bB, lB, resB));
			output.push_back(ResIM(ba, la, resa));
			output.push_back(ResIM(bb, lb, resb));
		}else if (B>=b && b>=A && A>=a){
			output.push_back(ResIM(bA, lA, resA));
			output.push_back(ResIM(ba, la, resa));
		}else if (B>=A && A>=b && b>=a){
			output.push_back(ResIM(bA, lA, resA));
			output.push_back(ResIM(ba, la, resa));
			output.push_back(ResIM(bb, lb, resb));
		}else{//B>=A && A>=a && a>=b
			output.push_back(ResIM(bA, lA, resA));
			output.push_back(ResIM(ba, la, resa));
			output.push_back(ResIM(bb, lb, resb));
		}

		return output;
	}

	std::vector< std::vector<double> > Calculator::GetAllBlInvariantMasses(){
		std::vector< std::vector<double> > output;
		output.clear();
		std::vector<double> outputa, outputb;
		outputa.clear();
		outputb.clear();

		outputa.push_back(W1a().M());
		outputa.push_back(W2a().M());
		outputb.push_back(W1b().M());
		outputb.push_back(W2b().M());

		output.push_back(outputa);
		output.push_back(outputb);

		return output;
	}

	TLorentzVector Calculator::CalcPerp(const TLorentzVector &p4Vis,
			const long double phi){
		const long double x=cos(phi), y=sin(phi);
		const long double scale=(p4Vis.Px()*x+p4Vis.Py()*y)/(x*x+y*y);
		TLorentzVector perpOutput;
		perpOutput.SetXYZM(p4Vis.Px()-scale*x, p4Vis.Py()-scale*y, 0.0, p4Vis.M());
		return perpOutput;		
	}

	TLorentzVector Calculator::CalcPerp(const TLorentzVector &p4Vis,
			const TLorentzVector &p4Upstream_temp){
      TLorentzVector p4Upstream = p4Upstream_temp;
		if (p4Upstream.Pt()==0.0){
			p4Upstream.SetXYZM(1.0, 0.0, p4Upstream.Pz(), p4Upstream.M());
		}else{
			p4Upstream.SetXYZM(p4Upstream.Px()/p4Upstream.Pt(),
					p4Upstream.Py()/p4Upstream.Pt(),
					p4Upstream.Pz()/p4Upstream.Pt(),
					p4Upstream.M());
		}
		const long double scale=(p4Vis.Px()*p4Upstream.Px()+p4Vis.Py()*p4Upstream.Py())/(p4Upstream.Perp2());
		TLorentzVector perpOutput;
		perpOutput.SetXYZM(p4Vis.Px()-scale*(p4Upstream.Px()),
				p4Vis.Py()-scale*(p4Upstream.Py()), 0.0, p4Vis.M());
		return perpOutput;		
	}

	double Calculator::CalcMt2Perp(const TLorentzVector &vis1, const TLorentzVector &vis2,
			const TLorentzVector &child, const long double mc){
		const TLorentzVector upstream=-(vis1+vis2+child);
		const TLorentzVector perpVis1=CalcPerp(vis1, upstream);
		const TLorentzVector perpVis2=CalcPerp(vis2, upstream);
		const TLorentzVector perpChild=CalcPerp(child, upstream);
		return CalcMt2Analytically(perpVis1, perpVis2, perpChild, mc);
	}

	double Calculator::CalcMt2Perp(const TLorentzVector &vis1, const TLorentzVector &vis2,
			const TLorentzVector &child, const long double mc, const long double phi){
		const TLorentzVector perpVis1=CalcPerp(vis1, phi);
		const TLorentzVector perpVis2=CalcPerp(vis2, phi);
		const TLorentzVector perpChild=CalcPerp(child, phi);
		return CalcMt2Analytically(perpVis1, perpVis2, perpChild, mc);
	}

	double Calculator::CalcMct(const TLorentzVector &p4Vis1, const TLorentzVector &p4Vis2){
		const long double Et1 = PosSqrt(p4Vis1.Perp2()+p4Vis1.M2());
		const long double Et2 = PosSqrt(p4Vis2.Perp2()+p4Vis2.M2());

		return (double)(PosSqrt(p4Vis1.M2()+p4Vis2.M2()+2.0L*(Et1*Et2+p4Vis1.Px()*p4Vis2.Px()+p4Vis1.Py()*p4Vis2.Py())));
	}

	double Calculator::CalcMctPerp(const TLorentzVector &upstream,
			const TLorentzVector &p4Vis1in, const TLorentzVector &p4Vis2in){
		const TLorentzVector p4Vis1=CalcPerp(p4Vis1in, upstream);
		const TLorentzVector p4Vis2=CalcPerp(p4Vis2in, upstream);

		return CalcMct(p4Vis1, p4Vis2);
	}

	double Calculator::CalcMctPerp(const TLorentzVector &p4Vis1in,
			const TLorentzVector &p4Vis2in, const long double phi){
		const TLorentzVector p4Vis1=CalcPerp(p4Vis1in, phi);
		const TLorentzVector p4Vis2=CalcPerp(p4Vis2in, phi);

		return CalcMct(p4Vis1, p4Vis2);
	}

	long double Calculator::PosSqrt(const long double x){
		return sqrt(std::max(x,0.0L));
	}

	double Calculator::PosSqrt(const double x){
		return sqrt(std::max(x,0.0));
	}

	double Calculator::CalcMt2Numerically(const TLorentzVector &vis1, const TLorentzVector &vis2,
			const TLorentzVector &child, const long double mc){
		long double ma, pxa, mb, pxb, pyb, pxc, pyc;
		if (vis1.M()<vis2.M()){
			pxa=PosSqrt(vis1.Px()*vis1.Px()+vis1.Py()*vis1.Py());
			ma=vis1.M();
			mb=vis2.M();

			const long double cospart=vis1.Px()/pxa;
			const long double sinpart=vis1.Py()/pxa;

			pxb=vis2.Px()*cospart+vis2.Py()*sinpart;
			pyb=vis2.Py()*cospart-vis2.Px()*sinpart;

			pxc=child.Px()*cospart+child.Py()*sinpart;
			pyc=child.Py()*cospart-child.Px()*sinpart;
		}else{
			pxa=PosSqrt(vis2.Px()*vis2.Px()+vis2.Py()*vis2.Py());
			ma=vis2.M();
			mb=vis1.M();

			const long double cospart=vis2.Px()/pxa;
			const long double sinpart=vis2.Py()/pxa;

			pxb=vis1.Px()*cospart+vis1.Py()*sinpart;
			pyb=vis1.Py()*cospart-vis1.Px()*sinpart;

			pxc=child.Px()*cospart+child.Py()*sinpart;
			pyc=child.Py()*cospart-child.Px()*sinpart;
		}

		long double mt2min=mb+mc;
		long double mt2max=0.0L;

		//Find some bounds
		if (mb>0.0){
			const long double testx=pxc-pxb*mc/mb;
			const long double testy=pyc-pyb*mc/mb;

			const long double testmt2=PosSqrt(ma*ma+mc*mc+2.0L*(PosSqrt((ma*ma+pxa*pxa)*(mc*mc+testx*testx+testy*testy))-pxa*testx));

			const long double altmt2a=PosSqrt(ma*ma+mc*mc+2.0L*(PosSqrt((ma*ma+pxa*pxa)*(mc*mc+(pxc*pxc+pyc*pyc)*0.25L))-pxa*pxc*0.5L));
			const long double altmt2b=PosSqrt(mb*mb+mc*mc+2.0L*(PosSqrt((mb*mb+pxb*pxb+pyb*pyb)*(mc*mc+(pxc*pxc+pyc*pyc)*0.25L))-0.5L*(pxb*pxc+pyb*pyc)));

			mt2max=std::min(testmt2,std::max(altmt2a,altmt2b));
		}else{
			const long double altmt2a=PosSqrt(ma*ma+mc*mc+2.0L*(PosSqrt((ma*ma+pxa*pxa)*(mc*mc+(pxc*pxc+pyc*pyc)*0.25L))-pxa*pxc*0.5L));
			const long double altmt2b=PosSqrt(mb*mb+mc*mc+2.0L*(PosSqrt((mb*mb+pxb*pxb+pyb*pyb)*(mc*mc+(pxc*pxc+pyc*pyc)*0.25L))-0.5L*(pxb*pxc+pyb*pyc)));

			mt2max=std::max(altmt2a,altmt2b);
		}

		//Check balanced/unbalanced case
		if (mt2max<mt2min){
			return mt2min;
		}

		//Have to resort to numerical methods. Precompute some coefficients...
		long double mt2=mt2min+0.5L*(mt2max-mt2min);

		const long double A1=4.0L*ma*ma;
		const long double C1=4.0L*(ma*ma+pxa*pxa);

		const long double A2=4.0L*(mb*mb+pyb*pyb);
		const long double B2=-8.0L*pxb*pyb;
		const long double C2=4.0L*(mb*mb+pxb*pxb);

		const long double G1=-4.0L*A1*C1;
		const long double G2=B2*B2-4.0L*A2*C2;

		//Do the minimization
		while (mt2min<mt2 && mt2<mt2max){
			//Get remaining quadratic coefficients
			const long double D1=4.0L*pxa*(mc*mc+ma*ma-mt2*mt2);
			const long double F1=(2.0L*mc*pxa-mc*mc+ma*ma)*(2.0L*mc*pxa+mc*mc-ma*ma)+2.0L*(mc*mc+ma*ma)*mt2*mt2-mt2*mt2*mt2*mt2;

			const long double D2=4.0L*(2.0L*pxb*pyb*pyc-2.0L*pxc*pyb*pyb-2.0L*mb*mb*pxc-mc*mc*pxb-mb*mb*pxb)+4.0L*pxb*mt2*mt2;
			const long double E2=4.0L*(2.0L*pyb*pxb*pxc-2.0L*pyc*pxb*pxb-2.0L*mb*mb*pyc-mc*mc*pyb-mb*mb*pyb)+4.0L*pyb*mt2*mt2;
			const long double F2=4.0L*pxb*pxb*pyc*pyc+4.0L*mb*mb*pyc*pyc-8.0L*pxb*pxc*pyb*pyc+4.0L*mc*mc*pyb*pyc+4.0L*mb*mb*pyb*pyc+4.0L*pxc*pxc*pyb*pyb+4.0L*mc*mc*pyb*pyb+4.0L*mb*mb*pxc*pxc+4.0L*mc*mc*pxb*pxc+4.0L*mb*mb*pxb*pxc+4.0L*mc*mc*pxb*pxb-mc*mc*mc*mc+2.0L*mb*mb*mc*mc-mb*mb*mb*mb-2.0L*(2.0L*pyb*pyc+2.0L*pxb*pxc-mc*mc-mb*mb)*mt2*mt2-mt2*mt2*mt2*mt2;

			//Some intermediate variables...
			const long double H1=-4.0L*C1*D1;
			const long double H2=2.0L*(B2*E2-2.0L*C2*D2);

			const long double I1=-4.0L*C1*F1;
			const long double I2=E2*E2-4.0L*C2*F2;

			const long double J=C2*C2*G1+C1*C1*G2-C1*C1*B2*B2;
			const long double K=C2*C2*H1+C1*C1*H2-2.0L*C1*C1*B2*E2;
			const long double L=C2*C2*I1+C1*C1*I2-C1*C1*E2*E2;

			//Calculate Sturm sequence. Possible instability in projected cases
			const long double p04=J*J-4.0L*C1*C1*C2*C2*G1*G2;
			const long double p03=2.0L*J*K-4.0L*C1*C1*C2*C2*(G1*H2+G2*H1);
			const long double p02=2.0L*J*L+K*K-4.0L*C1*C1*C2*C2*(G1*I2+G2*I1+H1*H2);
			const long double p01=2.0L*K*L-4.0L*C1*C1*C2*C2*(H1*I2+H2*I1);
			const long double p00=L*L-4.0L*C1*C1*C2*C2*I1*I2;

			const long double p13=4.0L*p04;
			const long double p12=3.0L*p03;
			const long double p11=2.0L*p02;
			const long double p10=p01;

			short negroots=0;
			short posroots=0;

		        long double p22, p21, p20, p31, p30, p40;

			if (fabs(p04)>1e-7){//4th degree
				p22=(3.0L*p03*p03-8.0L*p02*p04)/(16.0L*p04);
				p21=(p02*p03-6.0L*p01*p04)/(8.0L*p04);
				p20=(p01*p03-16.0L*p00*p04)/(16.0L*p04);

				if (fabs(p22)>1e-7){
					p31=(p12*p21*p22+p13*p20*p22-p11*p22*p22-p13*p21*p21)/(p22*p22);
					p30=(p12*p20*p22-p10*p22*p22-p13*p20*p21)/(p22*p22);

					if (fabs(p31)>1e-7){
						p40=(p21*p30*p31-p20*p31*p31-p22*p30*p30)/(p31*p31);
					}else{
						p31=0.0;
						p40=0.0;
					}
				}else{
					p22=0.0;
					p31=0.0;
					p40=0.0;
					if (fabs(p21)>1e-7){
						p30=(p11*p20*p21*p21+p13*p20*p20*p20-p10*p21*p21*p21-p12*p20*p20*p21)/(p21*p21*p21);
					}else{
						p21=0.0;
						p30=0.0;
					}
				}
					
				if ((p04<0.0 && p13<0.0) || (p04>0.0 && p13>0.0)){
					negroots+=1;
				}else if ((p04<0.0 && p13>0.0) || (p04>0.0 && p13<0.0)){
					posroots+=1;
				}
				
				if ((p13<0.0 && p22<0.0) || (p13>0.0 && p22>0.0)){
					negroots+=1;
				}else if ((p13<0.0 && p22>0.0) || (p13>0.0 && p22<0.0)){
					posroots+=1;
				}

				if ((p22<0.0 && p31<0.0) || (p22>0.0 && p31>0.0)){
					negroots+=1;
				}else if ((p22<0.0 && p31>0.0) || (p22>0.0 && p31<0.0)){
					posroots+=1;
				}

				if ((p31<0.0 && p40<0.0) || (p31>0.0 && p40>0.0)){
					negroots+=1;
				}else if ((p31<0.0 && p40>0.0) || (p31>0.0 && p40<0.0)){
					posroots+=1;
				}
				
			}else{
				if (fabs(p03)>1e-7){//3rd degree, not sure that this can happen...
					p21=2.0L*(p02*p02-3.0L*p01*p03)/(9.0L*p03);
					p20=(p01*p02-9.0L*p00*p03)/(9.0L*p03);
					if (fabs(p21)>1e-7){
						p30=(p11*p20*p21-p10*p21*p21-p12*p20*p20)/(p21*p21);
					}else{
						p21=0.0;
						p30=0.0;
					
						if ((p03<0.0 && p12<0.0) || (p03>0.0 && p12>0.0)){
							negroots+=1;
						}else if ((p03<0.0 && p12>0.0) || (p03>0.0 && p12<0.0)){
							posroots+=1;
						}
						
						if ((p12<0.0 && p21<0.0) || (p12>0.0 && p21>0.0)){
							negroots+=1;
						}else if ((p12<0.0 && p21>0.0) || (p12>0.0 && p21<0.0)){
							posroots+=1;
						}
							
						if ((p21<0.0 && p30<0.0) || (p21>0.0 && p30>0.0)){
							negroots+=1;
						}else if ((p21<0.0 && p30>0.0) || (p21>0.0 && p30<0.0)){
							posroots+=1;
						}
					}
				}else{
					if (fabs(p02)>1e-7){//2nd degree, perp/par go here?
						p20=(p01*p01-4.0L*p00*p02)/(4.0L*p02);

						if ((p02<0.0 && p11<0.0) || (p02>0.0 && p11>0.0)){
							negroots+=1;
						}else if ((p02<0.0 && p11>0.0) || (p02>0.0 && p11<0.0)){
							posroots+=1;
						}
						
						if ((p11<0.0 && p20<0.0) || (p11>0.0 && p20>0.0)){
							negroots+=1;
						}else if ((p11<0.0 && p20>0.0) || (p11>0.0 && p20<0.0)){
							posroots+=1;
						}
						
					}else{
						return mt2;//Should probably check these cases
					}
				}
			}
				
			//Number of solutions is difference between posroots && negroots
			if (posroots==negroots){
				mt2min=mt2;
			}else{
				mt2max=mt2;
			}

			mt2=mt2min+0.5L*(mt2max-mt2min);
		}								
		return (double)(mt2);
	}

	double Calculator::CalcMt2Analytically(const TLorentzVector &vis1,
			const TLorentzVector &vis2, const TLorentzVector &child,
			const long double mc){
		long double ma, pxa, mb, pxb, pyb, pxc, pyc;
		if (vis1.M()<vis2.M()){
			pxa=PosSqrt(vis1.Px()*vis1.Px()+vis1.Py()*vis1.Py());
			ma=vis1.M();
			mb=vis2.M();
		
			const long double cospart=vis1.Px()/pxa;
			const long double sinpart=vis1.Py()/pxa;
		
			pxb=vis2.Px()*cospart+vis2.Py()*sinpart;
			pyb=vis2.Py()*cospart-vis2.Px()*sinpart;
		
			pxc=child.Px()*cospart+child.Py()*sinpart;
			pyc=child.Py()*cospart-child.Px()*sinpart;		
		}else{
			pxa=PosSqrt(vis2.Px()*vis2.Px()+vis2.Py()*vis2.Py());
			ma=vis2.M();
			mb=vis1.M();
		
			const long double cospart=vis2.Px()/pxa;
			const long double sinpart=vis2.Py()/pxa;
		
			pxb=vis1.Px()*cospart+vis1.Py()*sinpart;
			pyb=vis1.Py()*cospart-vis1.Px()*sinpart;
		
			pxc=child.Px()*cospart+child.Py()*sinpart;
			pyc=child.Py()*cospart-child.Px()*sinpart;
		}
	
		long double mt2min=mb+mc;
		long double mt2max=0.0L;
	
		//Find some bounds
		if (mb>0.0){
			const long double testx=pxc-pxb*mc/mb;
			const long double testy=pyc-pyb*mc/mb;
		
			mt2max=PosSqrt(ma*ma+mc*mc+2.0L*(PosSqrt((ma*ma+pxa*pxa)*(mc*mc+testx*testx+testy*testy))-pxa*testx));
		}else{
			const long double altmt2a=PosSqrt(ma*ma+mc*mc+2.0L*(PosSqrt((ma*ma+pxa*pxa)*(mc*mc+(pxc*pxc+pyc*pyc)*0.25L))-pxa*pxc*0.5L));
			const long double altmt2b=PosSqrt(mb*mb+mc*mc+2.0L*(PosSqrt((mb*mb+pxb*pxb+pyb*pyb)*(mc*mc+(pxc*pxc+pyc*pyc)*0.25L))-0.5L*(pxb*pxc+pyb*pyc)));
		
			mt2max=std::max(altmt2a,altmt2b);
		}
		
		//Check balaned/unbalanced case
		if (mt2max<mt2min){
			return mt2min;
		}

		const long double AT=PosSqrt((vis1.M2()+vis1.Px()*vis1.Px()+vis1.Py()*vis1.Py())*(vis2.M2()+vis2.Px()*vis2.Px()+vis2.Py()*vis2.Py()))+vis1.Px()*vis2.Px()+vis1.Py()*vis2.Py();
		return (double)(PosSqrt(mc*mc+AT+PosSqrt((1.0L+4.0L*mc*mc/(2.0L*AT-vis1.M2()-vis2.M2()))*(AT*AT-vis1.M2()*vis2.M2()))));
	}
}
