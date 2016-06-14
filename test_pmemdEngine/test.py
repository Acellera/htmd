import htmd
adapt = htmd.AdaptiveRun()
adapt.nmin = 2
adapt.nmax = 3
adapt.nepochs = 2
adapt.ticadim = 0
adapt.metricsel1 = 'name CA'
adapt.generatorspath = './'
adapt.app = htmd.apps.pmemdlocal.PmemdLocal()
adapt.run()
