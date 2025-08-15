import pathsim
from pathsim import Simulation, Connection
import numpy as np
import matplotlib.pyplot as plt
import pathview

from tritium_model import (
    neutron_rate,
    total_irradiation_time,
    k_wall,
    k_top,
    baby_model,
    measured_TBR,
)

# Create global variables
A_OV = (baby_model.A_wall.to("m^2")).magnitude
A_IV = (baby_model.A_top.to("m^2")).magnitude
baby_vol = (baby_model.volume.to("m^3")).magnitude
k_IV = k_top.to("m/s").magnitude
k_OV = k_wall.to("m/s").magnitude
neutron_rate = neutron_rate.to("n/s").magnitude
irradiation_time = total_irradiation_time.to("s").magnitude


tbr = measured_TBR.to("dimensionless").magnitude

IV_gas_residence_time = 0.1 * baby_vol / (A_IV * k_IV)
OV_gas_residence_time = 2.2 * baby_vol / (A_OV * k_OV)

baby_residence_time = baby_vol / (k_IV * A_IV + k_OV * A_OV)

collection_efficiency = 0.95
conversion_efficiency = 1

# Create blocks
blocks, events = [], []

x_1_4 = pathsim.blocks.amplifier.Amplifier(gain=-1)
blocks.append(x_1_4)

neutron_rate_6 = pathsim.blocks.adder.Adder()
blocks.append(neutron_rate_6)

tbr_7 = pathsim.blocks.amplifier.Amplifier(gain=tbr)
blocks.append(tbr_7)

baby_8 = pathview.custom_pathsim_blocks.Process(
    residence_time=baby_residence_time,
)
blocks.append(baby_8)

iv_vial_activity_21 = pathsim.blocks.scope.Scope(
    labels=[
        "IV bubbler (vial1)",
        "IV bubbler (vial2)",
        "IV bubbler (vial3)",
        "IV bubbler (vial4)",
    ]
)
blocks.append(iv_vial_activity_21)

soluble_vs_insoluble_1 = pathview.custom_pathsim_blocks.Splitter2(f1=0.01, f2=0.99)
blocks.append(soluble_vs_insoluble_1)

iv_bubbler_23 = pathview.custom_pathsim_blocks.Bubbler(
    conversion_efficiency=conversion_efficiency,
    vial_efficiency=collection_efficiency,
    replacement_times=np.array([0.4, 0.6, 1, 1.5, 2.5, 4]) * 24 * 3600,
)
events_iv_bubbler_23 = iv_bubbler_23.create_reset_events()
events += events_iv_bubbler_23
blocks.append(iv_bubbler_23)

iv_vs_ov_24 = pathview.custom_pathsim_blocks.Splitter2(
    f1=A_IV * k_IV / (A_IV * k_IV + A_OV * k_OV),
    f2=A_OV * k_OV / (A_IV * k_IV + A_OV * k_OV),
)
blocks.append(iv_vs_ov_24)

soluble_vs_insoluble_25 = pathview.custom_pathsim_blocks.Splitter2(f1=0.01, f2=0.99)
blocks.append(soluble_vs_insoluble_25)

ov_bubbler_26 = pathview.custom_pathsim_blocks.Bubbler(
    conversion_efficiency=conversion_efficiency,
    vial_efficiency=collection_efficiency,
    replacement_times=np.array([1, 2.5, 4]) * 24 * 3600,
)
events_ov_bubbler_26 = ov_bubbler_26.create_reset_events()
events += events_ov_bubbler_26
blocks.append(ov_bubbler_26)

environment_27 = pathview.custom_pathsim_blocks.Integrator()
blocks.append(environment_27)

ov_vial_activity_28 = pathsim.blocks.scope.Scope(
    labels=[
        "OV bubbler (vial1)",
        "OV bubbler (vial2)",
        "OV bubbler (vial3)",
        "OV bubbler (vial4)",
    ]
)
blocks.append(ov_vial_activity_28)

baby_inventory_30 = pathsim.blocks.scope.Scope(labels=["BABY (inv)"])
blocks.append(baby_inventory_30)

stepsource_31 = pathsim.blocks.sources.StepSource(
    amplitude=neutron_rate, tau=irradiation_time
)
blocks.append(stepsource_31)

stepsource_32 = pathsim.blocks.sources.StepSource(amplitude=neutron_rate, tau=0)
blocks.append(stepsource_32)

neutron_source_33 = pathsim.blocks.scope.Scope(labels=["Neutron rate"])
blocks.append(neutron_source_33)

cumulative_release_34 = pathsim.blocks.scope.Scope(labels=["IV", "OV"])
blocks.append(cumulative_release_34)

iv_35 = pathview.custom_pathsim_blocks.Integrator()
blocks.append(iv_35)

ov_36 = pathview.custom_pathsim_blocks.Integrator()
blocks.append(ov_36)

iv_gas_37 = pathview.custom_pathsim_blocks.Process(
    residence_time=IV_gas_residence_time,
)
blocks.append(iv_gas_37)

ov_gas_38 = pathview.custom_pathsim_blocks.Process(
    residence_time=OV_gas_residence_time,
)
blocks.append(ov_gas_38)


# Create events


# Create connections

connections = [
    Connection(x_1_4[0], neutron_rate_6[0]),
    Connection(neutron_rate_6[0], tbr_7[0]),
    Connection(tbr_7[0], baby_8[0]),
    Connection(soluble_vs_insoluble_1["source1"], iv_bubbler_23["sample_in_soluble"]),
    Connection(soluble_vs_insoluble_1["source2"], iv_bubbler_23["sample_in_insoluble"]),
    Connection(iv_bubbler_23["vial1"], iv_vial_activity_21[0]),
    Connection(iv_bubbler_23["vial2"], iv_vial_activity_21[1]),
    Connection(iv_bubbler_23["vial3"], iv_vial_activity_21[2]),
    Connection(iv_bubbler_23["vial4"], iv_vial_activity_21[3]),
    Connection(baby_8["mass_flow_rate"], iv_vs_ov_24[0]),
    Connection(soluble_vs_insoluble_25["source1"], ov_bubbler_26["sample_in_soluble"]),
    Connection(
        soluble_vs_insoluble_25["source2"], ov_bubbler_26["sample_in_insoluble"]
    ),
    Connection(iv_bubbler_23["sample_out"], environment_27[0]),
    Connection(ov_bubbler_26["sample_out"], environment_27[1]),
    Connection(ov_bubbler_26["vial1"], ov_vial_activity_28[0]),
    Connection(ov_bubbler_26["vial2"], ov_vial_activity_28[1]),
    Connection(ov_bubbler_26["vial3"], ov_vial_activity_28[2]),
    Connection(ov_bubbler_26["vial4"], ov_vial_activity_28[3]),
    Connection(baby_8["inv"], baby_inventory_30[0]),
    Connection(stepsource_31[0], x_1_4[0]),
    Connection(stepsource_32[0], neutron_rate_6[1]),
    Connection(neutron_rate_6[0], neutron_source_33[0]),
    Connection(iv_35[0], cumulative_release_34[0]),
    Connection(ov_36[0], cumulative_release_34[1]),
    Connection(iv_vs_ov_24["source1"], iv_gas_37[0]),
    Connection(iv_gas_37["mass_flow_rate"], iv_35[0]),
    Connection(iv_gas_37["mass_flow_rate"], soluble_vs_insoluble_1[0]),
    Connection(iv_vs_ov_24["source2"], ov_gas_38[0]),
    Connection(ov_gas_38["mass_flow_rate"], ov_36[0]),
    Connection(ov_gas_38["mass_flow_rate"], soluble_vs_insoluble_25[0]),
]

# Create simulation
my_simulation = Simulation(
    blocks,
    connections,
    events=events,
    Solver=pathsim.solvers.SSPRK22,
    dt=100,
    dt_max=1.0,
    dt_min=1e-6,
    iterations_max=100,
    log=True,
    tolerance_fpi=1e-6,
    **{"tolerance_lte_rel": 1e-4, "tolerance_lte_abs": 1e-9},
)

if __name__ == "__main__":
    my_simulation.run(8 * 24 * 3600)

    # Optional: Plotting results
    scopes = [block for block in blocks if isinstance(block, pathsim.blocks.Scope)]
    fig, axs = plt.subplots(
        nrows=len(scopes), sharex=True, figsize=(10, 5 * len(scopes))
    )
    for i, scope in enumerate(scopes):
        plt.sca(axs[i] if len(scopes) > 1 else axs)
        time, data = scope.read()
        # plot the recorded data
        for p, d in enumerate(data):
            lb = scope.labels[p] if p < len(scope.labels) else f"port {p}"
            plt.plot(time, d, label=lb)
        plt.legend()
    plt.xlabel("Time")
    plt.show()
