import gen3gas as G3
import scipy.optimize 
import multiprocessing
import random


#-------------------------------
def log_entry(x, z, iternum, label):
    return "{:20s}Iter: {:04d}".format(label, iternum) + ("{:>10s}"*(len(x)+1)).format(*["{:.3f}".format(v) for v in [z]+x])

#-------------------------------
def f_eval(x, data):

    # z = sum([v*v for v in x])

    # return z
    x_unscaled = [x[i]*data.x_initial[i] for i in range(len(x))]

    data.variables.cycle_design_power,\
        data.variables.solar_multiple,\
        data.variables.dni_design_point,\
        data.variables.receiver_height,\
        data.variables.riser_inner_diam, \
        data.variables.downcomer_inner_diam, \
        data.variables.hours_tes,\
        data.variables.dT_approach_charge_hx,\
        data.variables.dT_approach_disch_hx = x_unscaled
    
    data.exec()
    lcoe = data.get_result_value('LCOE (real)')
    if lcoe < 0.1 or lcoe > 35:
        lcoe = float('nan')

    if lcoe < data.z_best['z']:
        data.z_best['z'] = lcoe
        data.z_best['xk'] = [v for v in x_unscaled]
        data.z_best['iter'] = data.current_iteration
    
    # logline = "{:20s}Iter: {:04d}".format(data.casename, data.current_iteration) + ("{:>10s}"*(len(x)+1)).format(*["{:.3f}".format(v) for v in [lcoe]+x_unscaled])
    logline = log_entry(x_unscaled, lcoe, data.current_iteration, data.casename)
    data.optimization_log += "\n" + logline

    print(logline)

    data.current_iteration += 1
    return lcoe

def f_callback(xk): #, data):
    # print(".... Iteration complete")
    xyz=None
    return

def optimize(thread_id):

    casename = ["north-bucket", "surround-bucket", "north-skip", "surround-skip"][thread_id % 4]
    case = "{:03d}_{:s}".format(thread_id, casename)

    # print(">>>> " + case)

    g = G3.Gen3opt()

    g.settings.is_north = 'north' in case
    g.optimization_log = ""
    g.settings.lift_technology = 'bucket' if 'bucket' in case else 'skip'
    g.casename = case 
    # g.variables.cycle_design_power = 50

    g.current_iteration = 0
    
    x0 = [
        g.variables.cycle_design_power      *(0.3 + random.random()*1.0),
        g.variables.solar_multiple          *(0.6 + random.random()*0.6),
        g.variables.dni_design_point        *(0.6 + random.random()*0.6),
        g.variables.receiver_height         *(0.5 + random.random()*1.5),
        g.variables.riser_inner_diam        *(0.6 + random.random()*1.0), 
        g.variables.downcomer_inner_diam    *(0.6 + random.random()*1.0), 
        g.variables.hours_tes               *(0.5 + random.random()*1.2),
        g.variables.dT_approach_charge_hx   *(0.4 + random.random()*2.0),
        g.variables.dT_approach_disch_hx    *(0.4 + random.random()*2.0),
    ]

    g.x_initial = [v for v in x0]
    g.z_best = {'z':float('inf'), 'xk':[v for v in x0], 'iter':-1}

    x0 = [1. for v in x0]

    scipy.optimize.fmin(f_eval, x0, args = ((g,)), ftol=0.01, xtol=0.01, maxfun=250, callback=f_callback) #, callback=f_update)

    logline = log_entry(g.z_best['xk'], g.z_best['z'], g.z_best['iter'], "***Best point:")
    g.optimization_log += "\n\n" + logline

    fout = open('optimization-log-'+case+'.txt', 'w')
    fout.write(g.optimization_log)
    fout.close()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    nthreads = 12
    all_args = [[i] for i in range(400)]

    pool = multiprocessing.Pool(processes=nthreads)
    pool.starmap(optimize, all_args)
    
    
    # optimize(0)
