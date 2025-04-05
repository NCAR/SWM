import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

def get_min_runtime_vs_n_proc(df):
    # Return Values of n_proc, fastest runtime for each n_proc, and run_idx for the fastest runtime for each n_proc
    n_proc_values = sorted(df['n_proc'].unique().tolist())
    t_min = []
    run_idx = []
    
    for n_proc in n_proc_values:
        np_mask = (df['n_proc'] == n_proc)
    
        fastest_t_for_this_np = df[np_mask]['runtime_max'].min()
        run_idx_of_fastest_t_for_this_np = int(df[np_mask & (df['runtime_max'] == fastest_t_for_this_np)]['run_idx'].values[0])
    
        t_min.append(fastest_t_for_this_np)    
        run_idx.append(run_idx_of_fastest_t_for_this_np)
    
    return n_proc_values, t_min, run_idx

def get_max_runtime_vs_n_proc(df):
    # Return Values of n_proc, slowest runtime for each n_proc, and run_idx for the fastest runtime for each n_proc
    n_proc_values = sorted(df['n_proc'].unique().tolist())
    t_min = []
    run_idx = []
    
    for n_proc in n_proc_values:
        np_mask = (df['n_proc'] == n_proc)
    
        slowest_t_for_this_np = df[np_mask]['runtime_max'].max()
        run_idx_of_slowest_t_for_this_np = int(df[np_mask & (df['runtime_max'] == slowest_t_for_this_np)]['run_idx'].values[0])
    
        t_min.append(slowest_t_for_this_np)    
        run_idx.append(run_idx_of_slowest_t_for_this_np)
    
    return n_proc_values, t_min, run_idx

def speedup(t_serial, t_parallel):
    return t_serial/np.array(t_parallel)

def efficiency(t_serial, t_parallel, n_proc):
    return speedup(t_serial, t_parallel)/np.array(n_proc)


def plot_multiple_timers_runtime(strong_scaling_results_dir, n_proc_values, run_index_values, fig):
    file_basenames = ['exclusive_runtime_UpdateIntermediateVariables',
                      'exclusive_runtime_UpdateNewVariables',
                      'exclusive_runtime_UpdateOldVariables',
                      'inclusive_runtime_FabArray::FillBoundary'
                     ]
    files = []
    for file_basename in file_basenames:
        file_path = os.path.join(strong_scaling_results_dir, file_basename+'.csv')
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
        files.append(file_path)

    plt.figure(fig.number) 

    for file_path in files:

      column_names = ['n_proc', 'run_idx', 'runtime_min', 'runtime_avg', 'runtime_max', 'max_percent']
      df = pd.read_csv(file_path, sep=',', header=0, names=column_names, index_col=False, skiprows=0, engine='python')

      time = []
      
      for (n_proc_value, fastest_run_idx_value) in zip(n_proc_values, run_index_values):
          # Get the data for the fastest running case for this specific value of n_proc
          mask = (df['n_proc'] == n_proc_value) & (df['run_idx'] == fastest_run_idx_value)
        
          time.append(float(df[mask]['runtime_max'].values[0]))
        
      # Will be used as x axis
      plt.plot(n_proc_values, time, label=os.path.basename(file_path))

def main():
    ###############################################################################
    # User Input
    ###############################################################################
    strong_scaling_results_dir = '/home/lalo/SWM/swm_AMReX/scaling_output_test'
    
    # Controls how the plots looks
    label_font_size=18
    figure_size=(8,8)
    
    ############################################################################
    # 
    ############################################################################

    if not os.path.isdir(strong_scaling_results_dir):
        raise NotADirectoryError(f"Directory not found: {strong_scaling_results_dir}")
    
    column_names = ['n_proc', 'run_idx', 'runtime_min', 'runtime_avg', 'runtime_max']
    df = pd.read_csv(os.path.join(strong_scaling_results_dir, 'total_runtime.csv'), sep=',', header=0, names=column_names, index_col=False, skiprows=0, engine='python')
    
    #print(df.to_string())
    
    all_runtimes_figure = plt.figure(figsize=figure_size)
    runtime_figure = plt.figure(figsize=figure_size)
    runtime_multi_timer_figure = plt.figure(figsize=figure_size)
    speedup_figure = plt.figure(figsize=figure_size)
    efficiency_figure = plt.figure(figsize=figure_size)
    
    ############################################################################
    # 
    ############################################################################
    
    plt.figure(all_runtimes_figure.number)
    plt.xlabel(r'$N_{\mathrm{proc}}$',fontsize=label_font_size)
    plt.ylabel(r'Runtime [s]',fontsize=label_font_size)
    #plt.legend()
    plt.tight_layout()
    
    # Fastest run for each n_proc
    n_proc, t_min, fastest_run_idx = get_min_runtime_vs_n_proc(df)
    plt.plot(n_proc, t_min, label='Fastest Run')
    
    for (n_proc_value, t_min_value, fastest_run_idx_value) in zip(n_proc, t_min, fastest_run_idx):
      print(f"Fastest run for n_proc {n_proc_value} was {t_min_value} [s] at run number {fastest_run_idx_value}")
    
    # Slowest run for each n_proc
    n_proc, t_max, slowest_run_idx = get_max_runtime_vs_n_proc(df)
    plt.plot(n_proc, t_max, label='Slowest Run')
    
    # Fill between the fastest and slowest run for each n_proc
    plt.fill_between(n_proc, t_min, t_max, alpha=0.2)
    
    # Plot x for all runs
    plt.plot(df['n_proc'].to_list(), df['runtime_max'].tolist(), 'xk', label='All Runs')
    
    plt.legend()
    
    ############################################################################
    # Runtime for fastest run only for each n_proc
    ############################################################################
    
    plt.figure(runtime_figure.number)
    
    plt.figure(runtime_figure.number)
    plt.xlabel(r'$N_{\mathrm{proc}}$',fontsize=label_font_size)
    plt.ylabel(r'Runtime [s]',fontsize=label_font_size)
    plt.tight_layout()
    
    plt.plot(n_proc, t_min, '-o', label='Actual')

    ############################################################################
    # Runtime for fastest run only for each n_proc.... but show multiple timers
    ############################################################################
    plt.figure(runtime_multi_timer_figure.number)
    plt.xlabel(r'$N_{\mathrm{proc}}$',fontsize=label_font_size)
    plt.ylabel(r'Runtime [s]',fontsize=label_font_size)
    plt.tight_layout()

    # Fastest run for each n_proc
    plt.plot(n_proc, t_min, label='toal runtime')
    plot_multiple_timers_runtime(strong_scaling_results_dir, n_proc, fastest_run_idx, runtime_multi_timer_figure)

    # Slowest run for each n_proc
    #plt.plot(n_proc, t_max, label='toal runtime')
    #plot_multiple_timers_runtime(strong_scaling_results_dir, n_proc, slowest_run_idx, runtime_multi_timer_figure)

    plt.legend()

    ############################################################################
    # Speedup
    ############################################################################
    
    plt.figure(speedup_figure.number)
    plt.plot([1,np.max(n_proc)], [1,np.max(n_proc)], '--k', label='Ideal') # Ideal speedup for reference
    plt.xlabel(r'$N_{\mathrm{proc}}$',fontsize=label_font_size)
    plt.ylabel(r'Speedup $= \frac{t_{\mathrm{serial}}}{t_\mathrm{parallel}}$',fontsize=label_font_size)
    plt.legend()
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    
    plt.plot(n_proc, speedup(t_min[0], t_min), '-o', label='Actual')
    
    ############################################################################
    # Efficiency
    ############################################################################
    
    plt.figure(efficiency_figure.number)
    plt.plot([1,np.max(n_proc)], [1,1], '--k', label='Ideal') # Ideal efficiency for reference
    plt.xlabel(r'$N_{\mathrm{proc}}$',fontsize=label_font_size)
    plt.ylabel(r'Efficiency $= \frac{\mathrm{Speedup}}{N_{\mathrm{proc}}} =  \frac{t_{\mathrm{serial}}}{N_{\mathrm{proc}}t_{\mathrm{parallel}}}$',fontsize=label_font_size)
    plt.legend()
    plt.tight_layout()
    
    plt.plot(n_proc, efficiency(t_min[0], t_min, n_proc), '-o', label='Actual')
    
    ############################################################################
    # Cleanup 
    ############################################################################
    
    plt.show()

if __name__ == "__main__":
    main()