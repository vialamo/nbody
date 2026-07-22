import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# 1. Load the HDF5 Data
file_path = "../build/data/CloudyData_UVB_HM2012.h5"

with h5py.File(file_path, 'r') as f:
    # Read Axes (remembering the specific formatting we found earlier)
    rho = f['/CoolingRates/Metals/Cooling'].attrs['Parameter1'][:]  # log10(nH)
    z = f['/CoolingRates/Metals/Cooling'].attrs['Parameter2'][:]    # linear z
    T = f['/CoolingRates/Metals/Cooling'].attrs['Temperature'][:]   # linear T (Kelvin)
    
    # Read 3D Data Grids. Shape is (Density, Redshift, Temperature)
    cool_prim = f['/CoolingRates/Primordial/Cooling'][:]
    heat_prim = f['/CoolingRates/Primordial/Heating'][:]
    cool_metal = f['/CoolingRates/Metals/Cooling'][:]
    heat_metal = f['/CoolingRates/Metals/Heating'][:]

# 2. Setup the Plot
fig, ax = plt.subplots(figsize=(10, 7))
plt.subplots_adjust(bottom=0.3)  # Make room for the sliders below

# Set initial slider indices (middle of the arrays)
init_d = len(rho) // 2
init_t = len(T) // 2

# Plot the 4 initial lines
line_cp, = ax.plot(z, cool_prim[init_d, :, init_t], lw=2, label='Primordial Cooling', color='blue')
line_hp, = ax.plot(z, heat_prim[init_d, :, init_t], lw=2, label='Primordial Heating', color='orange')
line_cm, = ax.plot(z, cool_metal[init_d, :, init_t], lw=2, label='Metal Cooling', color='cyan', linestyle='--')
line_hm, = ax.plot(z, heat_metal[init_d, :, init_t], lw=2, label='Metal Heating', color='red', linestyle='--')

# Plot styling
ax.set_yscale('log') # Log scale because rates vary over many orders of magnitude
ax.set_xlabel('Redshift (z)')
ax.set_ylabel('Rate (erg / s / cm^3)')
ax.grid(True, which="both", ls="--", alpha=0.5)
ax.legend()

# 3. Setup the Sliders
ax_d = plt.axes([0.15, 0.15, 0.65, 0.03])
ax_t = plt.axes([0.15, 0.10, 0.65, 0.03])

# The sliders will output the integer index of the arrays
slider_d = Slider(ax_d, 'Density Index', 0, len(rho)-1, valinit=init_d, valstep=1)
slider_t = Slider(ax_t, 'Temp Index', 0, len(T)-1, valinit=init_t, valstep=1)

# 4. The Update Function (called every time a slider moves)
def update(val):
    d_idx = int(slider_d.val)
    t_idx = int(slider_t.val)
    
    # Update the Y data for all 4 lines
    line_cp.set_ydata(cool_prim[d_idx, :, t_idx])
    line_hp.set_ydata(heat_prim[d_idx, :, t_idx])
    line_cm.set_ydata(cool_metal[d_idx, :, t_idx])
    line_hm.set_ydata(heat_metal[d_idx, :, t_idx])
    
    # Dynamically update the title to show the actual physical values
    current_rho = rho[d_idx]
    current_T = T[t_idx]
    ax.set_title(f'Thermodynamic Rates vs Redshift\nlog10(nH) = {current_rho:.2f} | T = {current_T:.2e} K')
    
    # Auto-scale the Y-axis to fit the new data
    ax.relim()
    ax.autoscale_view(scalex=False, scaley=True)
    
    fig.canvas.draw_idle()

# Link the sliders to the update function
slider_d.on_changed(update)
slider_t.on_changed(update)

# Trigger the initial title generation
update(None)

plt.show()