import yt
#from yt.frontends import boxlib
#from yt.frontends.boxlib.data_structures import AMReXDataset
from yt.frontends import amrex
from yt.frontends.amrex.data_structures import AMReXDataset

#yt.toggle_interactivity()

#ds = AMReXDataset("plt00000")
#ds = AMReXDataset("plt03900")

for filename in ["plt00000", "plt03900"]:

  ds = AMReXDataset(filename)
  ds.field_list
  
  sl_phi_old = yt.AxisAlignedSlicePlot(ds, 'z', ('boxlib', 'psi'), origin="native", axes_unit="unitary")
  sl_p = yt.AxisAlignedSlicePlot(ds, 'z', ('boxlib', 'p'), origin="native", axes_unit="unitary")
  sl_u = yt.AxisAlignedSlicePlot(ds, 'z', ('boxlib', 'u'), origin="native", axes_unit="unitary")
  #sl_u = yt.AxisAlignedSlicePlot(ds, 'z', ('boxlib', 'u'), origin="native")
  sl_v = yt.AxisAlignedSlicePlot(ds, 'z', ('boxlib', 'v'), origin="native", axes_unit="unitary")
  
  #sl_u.set_cmap(field=("boxlib", "u"), cmap="inferno")
  
  #sl.plots[("boxlib", "output_array")].axes.grid(True)
  
  sl_phi_old.save()
  sl_p.save()
  sl_u.save()
  sl_v.save()