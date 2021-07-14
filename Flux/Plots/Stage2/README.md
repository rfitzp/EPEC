# Asymptote scripts (see http://asymptote.sourceforge.net)

## filename.asy
- Usage:	  asy filename
- Produces:   filename.eps
---

	psi.asy      ... plots psi (R, Z)
	psir.asy     ... plots interpolated dPsi/dR (R, Z)
	psiz.asy     ... plots interpolated dPsi/dZ (R, Z)
	psirr.asy    ... plots interpolated d^2Psi/dR^2 (R, Z)
	psirz.asy    ... plots interpolated d^2Psi/dRdZ (R, Z)	
	psizz.asy    ... plots interpolated d^2Psi/dZ^2 (R, Z)

	qr.asy       ... plots q(r) profile
	qredge.asy   ... plots q(r) profile close to plasma edge
	qedge.asy    ... plots q(r) profile close to plasma edge (no comparison w. gFile)
	pr.asy       ... plots p(r) profile
	g.asy        ... plots g(r) profile
	r.asy        ... plots R(r) profile on midplane
	bt.asy       ... plots B_toroidal(r) on midplane
	bp.asy       ... plots B_poloidal(r) on midplane
	
	qgp.asy      ... plots (q/g)(PsiN)
	qp.asy       ... plots q(PsiN)
	qp1.asy      ... plots 1st derivative of q(PsiN)
	qp2.asy      ... plots 2nd derivative of q(PsiN)
	a1.asy       ... plots A1(PsiN)
	a2.asy       ... plots A2(PsiN)

	rational.asy ... plots rational surfaces

	rst.asy      ... plots R versus theta on rational surfaces
	zst.asy      ... plots Z versus theta on rational surfaces

	rnc.asy      ... plots R versus Theta on rational surfaces
	znc.asy      ... plots Z versus Theta on rational surfaces
	bnc.asy      ... plots B versus Theta on rational surfaces
	cnc.asy      ... plots dBdTheta versus Theta on rational surfaces
