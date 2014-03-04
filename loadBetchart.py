#!/bin/env python

import sys
import ROOT as r
r.gROOT.SetBatch(1)
r.gROOT.LoadMacro('vecUtils.h'+'+')
lv = r.Math.LorentzVector(r.Math.PtEtaPhiE4D('float'))
#sys.path.append('./analytic-nu')
import nuSolutions as ns
doubleNeutrinoSolutionsCheckLinAlg = ns.doubleNeutrinoSolutionsCheckLinAlg

