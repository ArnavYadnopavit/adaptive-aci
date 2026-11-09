# adaptive-aci

This repo is for the course project of Communication Systems EE3700, Adaptive Communication under Adjacent-Channel Interference.

The file `Baseband_randomstr_with_aci.m` contains the code of ACI-type `randomstr` and `ratioed`.

To change the type of ACI, search for `interference_type`.

In Ratioed ACI-Type, ACI happens for 'x' continuous blocks and then for 'y' blocks ACI doesn't happen. We can select the number of blocks with ACI interference.

To change the values of 'x' and 'y', search for `aci_block_pattern` in line 24 of the code.

In Randomstr ACI-Type, ACI occurs randomly in bursts of consecutive blocks. We can configure the maximum ACI duty cycle and the maximum number of continuous blocks for which ACI can remain active.

To change these values, search for `randomstr_config` in the code. The first value specifies the maximum duty cycle (%) and the second value specifies the maximum continuous ACI-ON burst length.

For example, `[40 5]` means that ACI can be active for a maximum of 40% of the total blocks, with no more than 5 consecutive ACI-ON blocks.
