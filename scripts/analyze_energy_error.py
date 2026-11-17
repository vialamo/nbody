#!/usr/bin/env python3
import h5py
import numpy as np
import glob
import os
import sys
import argparse
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider

def load_run(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    if not files:
        print(f"No HDF5 files found in {snapshot_dir}")
        return None

    times = []
    ke, ie, te = [], [], []
    w_grav, w_exp, e_switch, e_floor = [], [], [], []

    for fname in files:
        with h5py.File(fname, 'r') as f:
            t = f['Header'].attrs['simulation_time']
            gas = f['Gas']
            
            mass = gas['mass'][:]
            vx = gas['velocity_x'][:]
            vy = gas['velocity_y'][:]
            vz = gas['velocity_z'][:]
            u = gas['internal_energy'][:]
            
            K = np.sum(0.5 * mass * (vx**2 + vy**2 + vz**2))
            U = np.sum(mass * u)
            
            times.append(t)
            ke.append(K)
            ie.append(U)
            te.append(K + U)
            
            w_grav.append(gas.attrs.get('cumulative_gravitational_work', 0.0))
            w_exp.append(gas.attrs.get('cumulative_expansion_work', 0.0))
            e_switch.append(gas.attrs.get('cumulative_entropy_switch_energy', 0.0))
            e_floor.append(gas.attrs.get('cumulative_floor_injected_energy', 0.0))

    data = {
        'time': np.array(times),
        'ke': np.array(ke),
        'ie': np.array(ie),
        'te': np.array(te),
        'w_grav': np.array(w_grav),
        'w_exp': np.array(w_exp),
        'e_switch': np.array(e_switch),
        'e_floor': np.array(e_floor)
    }
    
    # Calculate Adjusted Energy (TE minus external work/injections)
    data['e_adjusted'] = (data['te'] - data['w_grav'] + data['w_exp'] 
                          - data['e_switch'] - data['e_floor'])
    
    # Fractional Error
    E0 = data['e_adjusted'][0]
    data['e_error'] = (data['e_adjusted'] - E0) / E0 if E0 != 0 else np.zeros_like(times)
    
    # Drift components (change relative to t=0)
    data['delta_ke'] = data['ke'] - data['ke'][0]
    data['delta_ie'] = data['ie'] - data['ie'][0]
    
    # Calculate the sum for visualization
    data['delta_sum'] = data['delta_ke'] + data['delta_ie']
    
    return data

def plot_energy_diagnostics_interactive(gt_dir, ipt_dir):
    print(f"Loading Global Timestep (GT) run: {gt_dir}")
    gt = load_run(gt_dir)
    
    print(f"Loading Individual Particle Timestep (IPT) run: {ipt_dir}")
    ipt = load_run(ipt_dir)
    
    if gt is None or ipt is None:
        sys.exit(1)
        
    max_idx = min(len(gt['time']), len(ipt['time'])) - 1

    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1], figure=fig)
    plt.subplots_adjust(bottom=0.15, top=0.92, wspace=0.25, hspace=0.30)
    fig.canvas.manager.set_window_title("Sedov Blastwave Energy Diagnostics")
    
    ax_err = fig.add_subplot(gs[0, 0])
    ax_drift = fig.add_subplot(gs[0, 1])
    ax_inj = fig.add_subplot(gs[1, 0])
    ax_work = fig.add_subplot(gs[1, 1])

    # 1. Total Fractional Energy Error
    line_err_gt, = ax_err.plot([], [], label='Global Timestep (GT)', color='tab:blue', lw=2)
    line_err_ipt, = ax_err.plot([], [], label='Individual Timestep (IPT)', color='tab:red', lw=2, linestyle='--')
    ax_err.set_title("Fractional Total Energy Error")
    ax_err.set_xlabel("Time")
    ax_err.set_ylabel(r"$\Delta E / E_0$")
    ax_err.grid(True, linestyle=':', alpha=0.6)
    ax_err.legend()
    
    # 2. Kinetic vs Internal Energy Drift
    line_ke_ipt, = ax_drift.plot([], [], label=r'$\Delta$ Kinetic Energy', color='tab:green', lw=2)
    line_ie_ipt, = ax_drift.plot([], [], label=r'$\Delta$ Internal Energy', color='tab:orange', lw=2)
    line_sum_ipt, = ax_drift.plot([], [], label=r'$\Sigma$ (KE + IE) (IPT)', color='tab:red', lw=2, linestyle='-.')
    
    line_ke_gt, = ax_drift.plot([], [], label=r'GT $\Delta$ KE (Baseline)', color='tab:green', lw=1, alpha=0.5)
    line_ie_gt, = ax_drift.plot([], [], label=r'GT $\Delta$ IE (Baseline)', color='tab:orange', lw=1, alpha=0.5)
    line_sum_gt, = ax_drift.plot([], [], label=r'GT $\Sigma$ (Baseline)', color='tab:blue', lw=1, linestyle=':')
    
    ax_drift.set_title("Energy Component Drift (Relative to t=0)")
    ax_drift.set_xlabel("Time")
    ax_drift.set_ylabel("Energy [Code Units]")
    ax_drift.grid(True, linestyle=':', alpha=0.6)
    ax_drift.legend()
    
    # 3. Artificial Energy Injections
    line_switch_ipt, = ax_inj.plot([], [], label='Entropy Switch Injection (IPT)', color='tab:purple', lw=2)
    line_floor_ipt, = ax_inj.plot([], [], label='Temp Floor Injection (IPT)', color='tab:brown', lw=2)
    line_switch_gt, = ax_inj.plot([], [], label='Entropy Switch (GT)', color='tab:purple', lw=1, linestyle='--')
    line_floor_gt, = ax_inj.plot([], [], label='Temp Floor (GT)', color='tab:brown', lw=1, linestyle='--')
    ax_inj.set_title("Artificial Energy Injections (Dual Energy / Floor)")
    ax_inj.set_xlabel("Time")
    ax_inj.set_ylabel("Cumulative Energy [Code Units]")
    ax_inj.grid(True, linestyle=':', alpha=0.6)
    ax_inj.legend()
    
    # 4. Net Work
    line_grav_ipt, = ax_work.plot([], [], label='Gravitational Work (IPT)', color='gray', lw=2)
    line_exp_ipt, = ax_work.plot([], [], label='Expansion Work (IPT)', color='black', lw=2)
    line_grav_gt, = ax_work.plot([], [], label='Grav Work (GT)', color='gray', lw=1, linestyle='--')
    line_exp_gt, = ax_work.plot([], [], label='Exp Work (GT)', color='black', lw=1, linestyle='--')
    ax_work.set_title("External Work Terms")
    ax_work.set_xlabel("Time")
    ax_work.set_ylabel("Cumulative Work [Code Units]")
    ax_work.grid(True, linestyle=':', alpha=0.6)
    ax_work.legend()

    vl_err_gt = ax_err.axvline(0, color='tab:blue', linestyle='-', lw=1.5, alpha=0.5)
    vl_err_ipt = ax_err.axvline(0, color='tab:red', linestyle='--', lw=1.5, alpha=0.5)
    vl_drift_gt = ax_drift.axvline(0, color='tab:blue', linestyle='-', lw=1.5, alpha=0.5)
    vl_drift_ipt = ax_drift.axvline(0, color='tab:red', linestyle='--', lw=1.5, alpha=0.5)
    vl_inj_gt = ax_inj.axvline(0, color='tab:blue', linestyle='-', lw=1.5, alpha=0.5)
    vl_inj_ipt = ax_inj.axvline(0, color='tab:red', linestyle='--', lw=1.5, alpha=0.5)
    vl_work_gt = ax_work.axvline(0, color='tab:blue', linestyle='-', lw=1.5, alpha=0.5)
    vl_work_ipt = ax_work.axvline(0, color='tab:red', linestyle='--', lw=1.5, alpha=0.5)

    def update(val):
        idx = int(time_slider.val)
        
        idx_gt = idx if idx < len(gt['time']) else len(gt['time']) - 1
        idx_ipt = idx if idx < len(ipt['time']) else len(ipt['time']) - 1
            
        time_gt = gt['time'][idx_gt]
        time_ipt = ipt['time'][idx_ipt]

        fig.suptitle(f"Sedov Blastwave Energy Diagnostics: GT (t={time_gt:.4f}) vs IPT (t={time_ipt:.4f})", fontsize=16)

        line_err_gt.set_data(gt['time'][:idx_gt+1], gt['e_error'][:idx_gt+1])
        line_err_ipt.set_data(ipt['time'][:idx_ipt+1], ipt['e_error'][:idx_ipt+1])
        
        line_ke_ipt.set_data(ipt['time'][:idx_ipt+1], ipt['delta_ke'][:idx_ipt+1])
        line_ie_ipt.set_data(ipt['time'][:idx_ipt+1], ipt['delta_ie'][:idx_ipt+1])
        line_sum_ipt.set_data(ipt['time'][:idx_ipt+1], ipt['delta_sum'][:idx_ipt+1])
        
        line_ke_gt.set_data(gt['time'][:idx_gt+1], gt['delta_ke'][:idx_gt+1])
        line_ie_gt.set_data(gt['time'][:idx_gt+1], gt['delta_ie'][:idx_gt+1])
        line_sum_gt.set_data(gt['time'][:idx_gt+1], gt['delta_sum'][:idx_gt+1])
        
        line_switch_ipt.set_data(ipt['time'][:idx_ipt+1], ipt['e_switch'][:idx_ipt+1])
        line_floor_ipt.set_data(ipt['time'][:idx_ipt+1], ipt['e_floor'][:idx_ipt+1])
        line_switch_gt.set_data(gt['time'][:idx_gt+1], gt['e_switch'][:idx_gt+1])
        line_floor_gt.set_data(gt['time'][:idx_gt+1], gt['e_floor'][:idx_gt+1])
        
        line_grav_ipt.set_data(ipt['time'][:idx_ipt+1], ipt['w_grav'][:idx_ipt+1])
        line_exp_ipt.set_data(ipt['time'][:idx_ipt+1], ipt['w_exp'][:idx_ipt+1])
        line_grav_gt.set_data(gt['time'][:idx_gt+1], gt['w_grav'][:idx_gt+1])
        line_exp_gt.set_data(gt['time'][:idx_gt+1], gt['w_exp'][:idx_gt+1])

        vl_err_gt.set_xdata([time_gt, time_gt])
        vl_err_ipt.set_xdata([time_ipt, time_ipt])
        vl_drift_gt.set_xdata([time_gt, time_gt])
        vl_drift_ipt.set_xdata([time_ipt, time_ipt])
        vl_inj_gt.set_xdata([time_gt, time_gt])
        vl_inj_ipt.set_xdata([time_ipt, time_ipt])
        vl_work_gt.set_xdata([time_gt, time_gt])
        vl_work_ipt.set_xdata([time_ipt, time_ipt])
        
        max_time = max(np.nanmax(gt['time']), np.nanmax(ipt['time']))
        
        ax_err.set_xlim(0, max_time)
        all_err = np.concatenate([gt['e_error'][:idx_gt+1], ipt['e_error'][:idx_ipt+1]])
        if len(all_err) > 0:
            err_pad = max(abs(np.nanmax(all_err)), abs(np.nanmin(all_err))) * 0.1
            err_pad = err_pad if err_pad > 0 else 0.01
            ax_err.set_ylim(np.nanmin(all_err) - err_pad, np.nanmax(all_err) + err_pad)

        ax_drift.set_xlim(0, max_time)
        all_drift = np.concatenate([ipt['delta_ke'][:idx_ipt+1], ipt['delta_ie'][:idx_ipt+1], ipt['delta_sum'][:idx_ipt+1],
                                    gt['delta_ke'][:idx_gt+1], gt['delta_ie'][:idx_gt+1], gt['delta_sum'][:idx_gt+1]])
        if len(all_drift) > 0:
            drift_pad = (np.nanmax(all_drift) - np.nanmin(all_drift)) * 0.1
            drift_pad = drift_pad if drift_pad > 0 else 0.1
            ax_drift.set_ylim(np.nanmin(all_drift) - drift_pad, np.nanmax(all_drift) + drift_pad)
            
        ax_inj.set_xlim(0, max_time)
        all_inj = np.concatenate([ipt['e_switch'][:idx_ipt+1], ipt['e_floor'][:idx_ipt+1],
                                  gt['e_switch'][:idx_gt+1], gt['e_floor'][:idx_gt+1]])
        if len(all_inj) > 0:
            inj_pad = (np.nanmax(all_inj) - np.nanmin(all_inj)) * 0.1
            inj_pad = inj_pad if inj_pad > 0 else 0.1
            ax_inj.set_ylim(np.nanmin(all_inj) - inj_pad, np.nanmax(all_inj) + inj_pad)

        ax_work.set_xlim(0, max_time)
        all_work = np.concatenate([ipt['w_grav'][:idx_ipt+1], ipt['w_exp'][:idx_ipt+1],
                                   gt['w_grav'][:idx_gt+1], gt['w_exp'][:idx_gt+1]])
        if len(all_work) > 0:
            work_pad = (np.nanmax(all_work) - np.nanmin(all_work)) * 0.1
            work_pad = work_pad if work_pad > 0 else 0.1
            ax_work.set_ylim(np.nanmin(all_work) - work_pad, np.nanmax(all_work) + work_pad)

        fig.canvas.draw_idle()

    ax_slider = fig.add_axes([0.15, 0.05, 0.7, 0.03])
    time_slider = Slider(
        ax=ax_slider,
        label='Time Step Index',
        valmin=0,
        valmax=max_idx,
        valinit=max_idx,
        valfmt='%0.0f'
    )
    
    time_slider.on_changed(update)
    update(max_idx)
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare Sedov Energy Diagnostics interactively between GT and IPT.")
    parser.add_argument("gt_dir", type=str, help="Path to Global Timestep snapshot directory.")
    parser.add_argument("ipt_dir", type=str, help="Path to Individual Particle Timestep snapshot directory.")
    
    args = parser.parse_args()
    plot_energy_diagnostics_interactive(args.gt_dir, args.ipt_dir)

