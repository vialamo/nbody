import pandas as pd
import numpy as np
import glob
import os
import sys
import argparse
import matplotlib.pyplot as plt

def generate_performance_dashboard(sim_folder):
    csv_path = os.path.join(sim_folder, 'diagnostics.csv')
    if not os.path.exists(csv_path):
        print(f"Error: Could not find {csv_path}")
        return

    # Load diagnostics data
    df = pd.read_csv(csv_path)
    if df.empty:
        print("Error: diagnostics.csv is empty.")
        return

    print("Data extraction complete. Generating Matplotlib dashboard...")

    # Plotting Grid Setup
    fig, axs = plt.subplots(2, 3, figsize=(18, 10))
    fig.canvas.manager.set_window_title(f"Performance Dashboard: {os.path.basename(os.path.normpath(sim_folder))}")
    
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # ---------------------------------------------------------
    # (0, 0) Performance Summary Info Card
    # ---------------------------------------------------------
    axs[0, 0].axis('off')
    total_cycles = df['cycle'].iloc[-1] if not df.empty else 0
    total_wall_time = df['cumulative_wall_time'].iloc[-1] if 'cumulative_wall_time' in df.columns else 0
    peak_mem = df['memory_peak'].max() if 'memory_peak' in df.columns else 0
    avg_step_time = df['wall_time_total'].mean() if 'wall_time_total' in df.columns else 0
    
    info_text = (
        f"PERFORMANCE SUMMARY\n\n"
        f"Run Directory:\n{os.path.basename(os.path.normpath(sim_folder))}\n\n"
        f"{'Total Cycles':<18}: {total_cycles}\n"
        f"{'Total Wall Time':<18}: {total_wall_time:.2f} s\n"
        f"{'Avg Time/Step':<18}: {avg_step_time:.4f} s\n"
        f"{'Peak Memory':<18}: {peak_mem:.2f} MB"
    )
    
    axs[0, 0].text(0.5, 0.5, info_text, fontsize=12, color='black', va='center', ha='center',
                   multialignment='left', bbox=dict(boxstyle='round,pad=1.0', facecolor='#f8f9fa', alpha=1.0, edgecolor='#cccccc'), family='monospace')

    # ---------------------------------------------------------
    # (0, 1) Execution Time Breakdown (Stacked Area Chart)
    # ---------------------------------------------------------
    time_cols = [
        'wall_time_tree', 'wall_time_density', 'wall_time_hydro', 
        'wall_time_pp', 'wall_time_pm', 'wall_time_cooling', 'wall_time_io'
    ]
    time_cols = [c for c in time_cols if c in df.columns]
    
    if time_cols:
        labels = [c.replace('wall_time_', '').upper() for c in time_cols]
        axs[0, 1].stackplot(df['cycle'], *[df[c] for c in time_cols], labels=labels, colors=colors[:len(time_cols)], alpha=0.8)
        axs[0, 1].set(xlabel='Cycle', ylabel='Wall Time (s)', title='Execution Time Breakdown')
        axs[0, 1].legend(loc='upper left', fontsize=8)
        axs[0, 1].grid(True, linestyle='--', alpha=0.6)
    else:
        axs[0, 1].text(0.5, 0.5, 'No Timing Data', ha='center', va='center')

    # ---------------------------------------------------------
    # (0, 2) Timestep Constraints Evolution (Log Scale)
    # ---------------------------------------------------------
    dt_cols = ['dt_final', 'dt_cfl', 'dt_gravity', 'dt_cool']
    dt_cols = [c for c in dt_cols if c in df.columns]
    
    if dt_cols:
        for idx, col in enumerate(dt_cols):
            axs[0, 2].plot(df['sim_time'], df[col], label=col, color=colors[idx % len(colors)], lw=2)
        axs[0, 2].set(yscale='log', xlabel='Simulation Time', ylabel='Timestep Size (dt)', title='Timestep Constraints Evolution')
        axs[0, 2].legend(loc='upper right', fontsize=8)
        axs[0, 2].grid(True, which="both", ls="--", alpha=0.6)
    else:
        axs[0, 2].text(0.5, 0.5, 'No Timestep Data', ha='center', va='center')

    # ---------------------------------------------------------
    # (1, 0) Memory Usage Tracker
    # ---------------------------------------------------------
    mem_cols = ['memory_peak', 'memory_current']
    mem_cols = [c for c in mem_cols if c in df.columns]
    
    if mem_cols:
        for idx, col in enumerate(mem_cols):
            label = col.replace('_', ' ').title()
            axs[1, 0].plot(df['cycle'], df[col], label=label, color=colors[idx % len(colors)], lw=2)
        axs[1, 0].set(xlabel='Cycle', ylabel='Memory (MB)', title='Resident Set Size (RSS) Memory Tracking')
        axs[1, 0].legend(loc='lower right', fontsize=8)
        axs[1, 0].grid(True, linestyle='--', alpha=0.6)
    else:
        axs[1, 0].text(0.5, 0.5, 'No Memory Data', ha='center', va='center')

    # ---------------------------------------------------------
    # (1, 1) Cumulative Wall Time
    # ---------------------------------------------------------
    if 'cumulative_wall_time' in df.columns:
        axs[1, 1].plot(df['cycle'], df['cumulative_wall_time'], color='crimson', lw=2, label='Cumulative Time')
        axs[1, 1].fill_between(df['cycle'], df['cumulative_wall_time'], color='crimson', alpha=0.2)
        axs[1, 1].set(xlabel='Cycle', ylabel='Total Wall Time (s)', title='Cumulative Simulation Wall Time')
        axs[1, 1].legend(loc='upper left', fontsize=8)
        axs[1, 1].grid(True, linestyle='--', alpha=0.6)
    else:
        axs[1, 1].text(0.5, 0.5, 'No Cumulative Time Data', ha='center', va='center')

    # ---------------------------------------------------------
    # (1, 2) GPU Performance Breakdown
    # ---------------------------------------------------------
    gpu_cols = ['wall_time_GPU_transf', 'wall_time_GPU_back', 'wall_time_GPU_compute']
    gpu_cols = [c for c in gpu_cols if c in df.columns]
    
    if gpu_cols and df[gpu_cols].sum().sum() > 0:
        labels = ['Transf (H->D)', 'Transf (D->H)', 'Compute']
        axs[1, 2].stackplot(df['cycle'], *[df[c] for c in gpu_cols], labels=labels, colors=colors[:len(gpu_cols)], alpha=0.8)
        axs[1, 2].set(xlabel='Cycle', ylabel='Wall Time (s)', title='GPU Operations Overhead')
        axs[1, 2].legend(loc='upper left', fontsize=8)
        axs[1, 2].grid(True, linestyle='--', alpha=0.6)
    else:
        axs[1, 2].text(0.5, 0.5, 'No GPU Data', ha='center', va='center')

    # ---------------------------------------------------------
    # Finalize Layout & Show
    # ---------------------------------------------------------
    plt.tight_layout(rect=[0, 0.01, 1, 0.99], h_pad=2.2, w_pad=2.0)
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a matplotlib performance report from simulation diagnostics.")
    
    parser.add_argument(
        "path", 
        type=str, 
        help="Path to the simulation folder containing diagnostics.csv"
    )
    
    parser.add_argument(
        "-l", "--latest", 
        action="store_true", 
        help="Automatically find and load the most recent in directory inside the provided path"
    )
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    target_dir = args.path

    if args.latest:
        search_pattern = os.path.join(target_dir, "run_*")
        runs = sorted(glob.glob(search_pattern))
        
        if runs:
            target_dir = runs[-1]
            print(f"Auto-selected latest run: {target_dir}")
        else:
            print(f"Error: No run directories found in '{args.path}'")
            sys.exit(1)
            
    generate_performance_dashboard(target_dir)
