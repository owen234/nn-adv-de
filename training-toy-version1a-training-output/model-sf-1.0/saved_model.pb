Јт
лО
^
AssignVariableOp
resource
value"dtype"
dtypetype"
validate_shapebool( 
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
8
Const
output"dtype"
valuetensor"
dtypetype
.
Identity

input"T
output"T"	
Ttype
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	

MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool("
allow_missing_filesbool( 

NoOp
M
Pack
values"T*N
output"T"
Nint(0"	
Ttype"
axisint 
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetype
E
Relu
features"T
activations"T"
Ttype:
2	
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
?
Select
	condition

t"T
e"T
output"T"	
Ttype
H
ShardedFilename
basename	
shard

num_shards
filename
С
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring Ј
@
StaticRegexFullMatch	
input

output
"
patternstring
ї
StridedSlice

input"T
begin"Index
end"Index
strides"Index
output"T"	
Ttype"
Indextype:
2	"

begin_maskint "
end_maskint "
ellipsis_maskint "
new_axis_maskint "
shrink_axis_maskint 
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 

VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 "serve*2.11.02unknown8Џ
^
countVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namecount
W
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
dtype0
^
totalVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nametotal
W
total/Read/ReadVariableOpReadVariableOptotal*
_output_shapes
: *
dtype0
b
count_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_1
[
count_1/Read/ReadVariableOpReadVariableOpcount_1*
_output_shapes
: *
dtype0
b
total_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_1
[
total_1/Read/ReadVariableOpReadVariableOptotal_1*
_output_shapes
: *
dtype0
~
Adam/v/dense_8/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*$
shared_nameAdam/v/dense_8/bias
w
'Adam/v/dense_8/bias/Read/ReadVariableOpReadVariableOpAdam/v/dense_8/bias*
_output_shapes
:*
dtype0
~
Adam/m/dense_8/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*$
shared_nameAdam/m/dense_8/bias
w
'Adam/m/dense_8/bias/Read/ReadVariableOpReadVariableOpAdam/m/dense_8/bias*
_output_shapes
:*
dtype0

Adam/v/dense_8/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*&
shared_nameAdam/v/dense_8/kernel

)Adam/v/dense_8/kernel/Read/ReadVariableOpReadVariableOpAdam/v/dense_8/kernel*
_output_shapes

:2*
dtype0

Adam/m/dense_8/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*&
shared_nameAdam/m/dense_8/kernel

)Adam/m/dense_8/kernel/Read/ReadVariableOpReadVariableOpAdam/m/dense_8/kernel*
_output_shapes

:2*
dtype0
~
Adam/v/dense_7/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*$
shared_nameAdam/v/dense_7/bias
w
'Adam/v/dense_7/bias/Read/ReadVariableOpReadVariableOpAdam/v/dense_7/bias*
_output_shapes
:2*
dtype0
~
Adam/m/dense_7/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*$
shared_nameAdam/m/dense_7/bias
w
'Adam/m/dense_7/bias/Read/ReadVariableOpReadVariableOpAdam/m/dense_7/bias*
_output_shapes
:2*
dtype0

Adam/v/dense_7/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*&
shared_nameAdam/v/dense_7/kernel

)Adam/v/dense_7/kernel/Read/ReadVariableOpReadVariableOpAdam/v/dense_7/kernel*
_output_shapes

:22*
dtype0

Adam/m/dense_7/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*&
shared_nameAdam/m/dense_7/kernel

)Adam/m/dense_7/kernel/Read/ReadVariableOpReadVariableOpAdam/m/dense_7/kernel*
_output_shapes

:22*
dtype0
~
Adam/v/dense_6/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*$
shared_nameAdam/v/dense_6/bias
w
'Adam/v/dense_6/bias/Read/ReadVariableOpReadVariableOpAdam/v/dense_6/bias*
_output_shapes
:2*
dtype0
~
Adam/m/dense_6/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*$
shared_nameAdam/m/dense_6/bias
w
'Adam/m/dense_6/bias/Read/ReadVariableOpReadVariableOpAdam/m/dense_6/bias*
_output_shapes
:2*
dtype0

Adam/v/dense_6/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*&
shared_nameAdam/v/dense_6/kernel

)Adam/v/dense_6/kernel/Read/ReadVariableOpReadVariableOpAdam/v/dense_6/kernel*
_output_shapes

:2*
dtype0

Adam/m/dense_6/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*&
shared_nameAdam/m/dense_6/kernel

)Adam/m/dense_6/kernel/Read/ReadVariableOpReadVariableOpAdam/m/dense_6/kernel*
_output_shapes

:2*
dtype0
n
learning_rateVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namelearning_rate
g
!learning_rate/Read/ReadVariableOpReadVariableOplearning_rate*
_output_shapes
: *
dtype0
f
	iterationVarHandleOp*
_output_shapes
: *
dtype0	*
shape: *
shared_name	iteration
_
iteration/Read/ReadVariableOpReadVariableOp	iteration*
_output_shapes
: *
dtype0	
p
dense_8/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_namedense_8/bias
i
 dense_8/bias/Read/ReadVariableOpReadVariableOpdense_8/bias*
_output_shapes
:*
dtype0
x
dense_8/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*
shared_namedense_8/kernel
q
"dense_8/kernel/Read/ReadVariableOpReadVariableOpdense_8/kernel*
_output_shapes

:2*
dtype0
p
dense_7/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*
shared_namedense_7/bias
i
 dense_7/bias/Read/ReadVariableOpReadVariableOpdense_7/bias*
_output_shapes
:2*
dtype0
x
dense_7/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*
shared_namedense_7/kernel
q
"dense_7/kernel/Read/ReadVariableOpReadVariableOpdense_7/kernel*
_output_shapes

:22*
dtype0
p
dense_6/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*
shared_namedense_6/bias
i
 dense_6/bias/Read/ReadVariableOpReadVariableOpdense_6/bias*
_output_shapes
:2*
dtype0
x
dense_6/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*
shared_namedense_6/kernel
q
"dense_6/kernel/Read/ReadVariableOpReadVariableOpdense_6/kernel*
_output_shapes

:2*
dtype0
z
serving_default_input_3Placeholder*'
_output_shapes
:џџџџџџџџџ*
dtype0*
shape:џџџџџџџџџ
Ђ
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_3dense_6/kerneldense_6/biasdense_7/kerneldense_7/biasdense_8/kerneldense_8/bias*
Tin
	2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *.
f)R'
%__inference_signature_wrapper_1184256

NoOpNoOp
Ѓ2
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*о1
valueд1Bб1 BЪ1
Ї
layer-0
layer-1
layer-2
layer_with_weights-0
layer-3
layer-4
	variables
trainable_variables
regularization_losses
		keras_api

__call__
*&call_and_return_all_conditional_losses
_default_save_signature
	optimizer

signatures*
* 

	keras_api* 

	keras_api* 

layer_with_weights-0
layer-0
layer_with_weights-1
layer-1
layer_with_weights-2
layer-2
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses*

	keras_api* 
.
0
1
2
3
4
 5*
.
0
1
2
3
4
 5*
* 
А
!non_trainable_variables

"layers
#metrics
$layer_regularization_losses
%layer_metrics
	variables
trainable_variables
regularization_losses

__call__
_default_save_signature
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses*
6
&trace_0
'trace_1
(trace_2
)trace_3* 
6
*trace_0
+trace_1
,trace_2
-trace_3* 
* 

.
_variables
/_iterations
0_learning_rate
1_index_dict
2
_momentums
3_velocities
4_update_step_xla*

5serving_default* 
* 
* 
І
6	variables
7trainable_variables
8regularization_losses
9	keras_api
:__call__
*;&call_and_return_all_conditional_losses

kernel
bias*
І
<	variables
=trainable_variables
>regularization_losses
?	keras_api
@__call__
*A&call_and_return_all_conditional_losses

kernel
bias*
І
B	variables
Ctrainable_variables
Dregularization_losses
E	keras_api
F__call__
*G&call_and_return_all_conditional_losses

kernel
 bias*
.
0
1
2
3
4
 5*
.
0
1
2
3
4
 5*
* 

Hnon_trainable_variables

Ilayers
Jmetrics
Klayer_regularization_losses
Llayer_metrics
	variables
trainable_variables
regularization_losses
__call__
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses*
6
Mtrace_0
Ntrace_1
Otrace_2
Ptrace_3* 
6
Qtrace_0
Rtrace_1
Strace_2
Ttrace_3* 
* 
NH
VARIABLE_VALUEdense_6/kernel&variables/0/.ATTRIBUTES/VARIABLE_VALUE*
LF
VARIABLE_VALUEdense_6/bias&variables/1/.ATTRIBUTES/VARIABLE_VALUE*
NH
VARIABLE_VALUEdense_7/kernel&variables/2/.ATTRIBUTES/VARIABLE_VALUE*
LF
VARIABLE_VALUEdense_7/bias&variables/3/.ATTRIBUTES/VARIABLE_VALUE*
NH
VARIABLE_VALUEdense_8/kernel&variables/4/.ATTRIBUTES/VARIABLE_VALUE*
LF
VARIABLE_VALUEdense_8/bias&variables/5/.ATTRIBUTES/VARIABLE_VALUE*
* 
'
0
1
2
3
4*

U0
V1*
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
b
/0
W1
X2
Y3
Z4
[5
\6
]7
^8
_9
`10
a11
b12*
SM
VARIABLE_VALUE	iteration0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUElearning_rate3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
* 
.
W0
Y1
[2
]3
_4
a5*
.
X0
Z1
\2
^3
`4
b5*
P
ctrace_0
dtrace_1
etrace_2
ftrace_3
gtrace_4
htrace_5* 
* 

0
1*

0
1*
* 

inon_trainable_variables

jlayers
kmetrics
llayer_regularization_losses
mlayer_metrics
6	variables
7trainable_variables
8regularization_losses
:__call__
*;&call_and_return_all_conditional_losses
&;"call_and_return_conditional_losses*

ntrace_0* 

otrace_0* 

0
1*

0
1*
* 

pnon_trainable_variables

qlayers
rmetrics
slayer_regularization_losses
tlayer_metrics
<	variables
=trainable_variables
>regularization_losses
@__call__
*A&call_and_return_all_conditional_losses
&A"call_and_return_conditional_losses*

utrace_0* 

vtrace_0* 

0
 1*

0
 1*
* 

wnon_trainable_variables

xlayers
ymetrics
zlayer_regularization_losses
{layer_metrics
B	variables
Ctrainable_variables
Dregularization_losses
F__call__
*G&call_and_return_all_conditional_losses
&G"call_and_return_conditional_losses*

|trace_0* 

}trace_0* 
* 

0
1
2*
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
:
~	variables
	keras_api

total

count*
M
	variables
	keras_api

total

count

_fn_kwargs*
`Z
VARIABLE_VALUEAdam/m/dense_6/kernel1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUEAdam/v/dense_6/kernel1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEAdam/m/dense_6/bias1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEAdam/v/dense_6/bias1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUEAdam/m/dense_7/kernel1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUEAdam/v/dense_7/kernel1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEAdam/m/dense_7/bias1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEAdam/v/dense_7/bias1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUEAdam/m/dense_8/kernel1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/v/dense_8/kernel2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/m/dense_8/bias2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/v/dense_8/bias2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUE*
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 

0
1*

~	variables*
UO
VARIABLE_VALUEtotal_14keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEcount_14keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE*

0
1*

	variables*
SM
VARIABLE_VALUEtotal4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUE*
SM
VARIABLE_VALUEcount4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUE*
* 
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
Е	
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename"dense_6/kernel/Read/ReadVariableOp dense_6/bias/Read/ReadVariableOp"dense_7/kernel/Read/ReadVariableOp dense_7/bias/Read/ReadVariableOp"dense_8/kernel/Read/ReadVariableOp dense_8/bias/Read/ReadVariableOpiteration/Read/ReadVariableOp!learning_rate/Read/ReadVariableOp)Adam/m/dense_6/kernel/Read/ReadVariableOp)Adam/v/dense_6/kernel/Read/ReadVariableOp'Adam/m/dense_6/bias/Read/ReadVariableOp'Adam/v/dense_6/bias/Read/ReadVariableOp)Adam/m/dense_7/kernel/Read/ReadVariableOp)Adam/v/dense_7/kernel/Read/ReadVariableOp'Adam/m/dense_7/bias/Read/ReadVariableOp'Adam/v/dense_7/bias/Read/ReadVariableOp)Adam/m/dense_8/kernel/Read/ReadVariableOp)Adam/v/dense_8/kernel/Read/ReadVariableOp'Adam/m/dense_8/bias/Read/ReadVariableOp'Adam/v/dense_8/bias/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOpConst*%
Tin
2	*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2 *0J 8 *)
f$R"
 __inference__traced_save_1184650
а
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenamedense_6/kerneldense_6/biasdense_7/kerneldense_7/biasdense_8/kerneldense_8/bias	iterationlearning_rateAdam/m/dense_6/kernelAdam/v/dense_6/kernelAdam/m/dense_6/biasAdam/v/dense_6/biasAdam/m/dense_7/kernelAdam/v/dense_7/kernelAdam/m/dense_7/biasAdam/v/dense_7/biasAdam/m/dense_8/kernelAdam/v/dense_8/kernelAdam/m/dense_8/biasAdam/v/dense_8/biastotal_1count_1totalcount*$
Tin
2*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2 *0J 8 *,
f'R%
#__inference__traced_restore_1184732эЇ
њ

)__inference_model_2_layer_call_fn_1184169
input_3
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityЂStatefulPartitionedCall
StatefulPartitionedCallStatefulPartitionedCallinput_3unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_model_2_layer_call_and_return_conditional_losses_1184137s
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*+
_output_shapes
:џџџџџџџџџ`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:џџџџџџџџџ
!
_user_specified_name	input_3
љ

.__inference_sequential_2_layer_call_fn_1184448

inputs
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityЂStatefulPartitionedCall
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *R
fMRK
I__inference_sequential_2_layer_call_and_return_conditional_losses_1183963o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџ`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
Й
P
$__inference__update_step_xla_1184399
gradient
variable:22*
_XlaMustCompile(*(
_construction_contextkEagerRuntime*
_input_shapes
:22: *
	_noinline(:H D

_output_shapes

:22
"
_user_specified_name
gradient:($
"
_user_specified_name
variable
Ч

)__inference_dense_8_layer_call_fn_1184545

inputs
unknown:2
	unknown_0:
identityЂStatefulPartitionedCallо
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_dense_8_layer_call_and_return_conditional_losses_1183873o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџ`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ2: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ2
 
_user_specified_nameinputs


I__inference_sequential_2_layer_call_and_return_conditional_losses_1183963

inputs!
dense_6_1183947:2
dense_6_1183949:2!
dense_7_1183952:22
dense_7_1183954:2!
dense_8_1183957:2
dense_8_1183959:
identityЂdense_6/StatefulPartitionedCallЂdense_7/StatefulPartitionedCallЂdense_8/StatefulPartitionedCallє
dense_6/StatefulPartitionedCallStatefulPartitionedCallinputsdense_6_1183947dense_6_1183949*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_dense_6_layer_call_and_return_conditional_losses_1183840
dense_7/StatefulPartitionedCallStatefulPartitionedCall(dense_6/StatefulPartitionedCall:output:0dense_7_1183952dense_7_1183954*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_dense_7_layer_call_and_return_conditional_losses_1183857
dense_8/StatefulPartitionedCallStatefulPartitionedCall(dense_7/StatefulPartitionedCall:output:0dense_8_1183957dense_8_1183959*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_dense_8_layer_call_and_return_conditional_losses_1183873w
IdentityIdentity(dense_8/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџЌ
NoOpNoOp ^dense_6/StatefulPartitionedCall ^dense_7/StatefulPartitionedCall ^dense_8/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 2B
dense_6/StatefulPartitionedCalldense_6/StatefulPartitionedCall2B
dense_7/StatefulPartitionedCalldense_7/StatefulPartitionedCall2B
dense_8/StatefulPartitionedCalldense_8/StatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs


ѕ
D__inference_dense_7_layer_call_and_return_conditional_losses_1184536

inputs0
matmul_readvariableop_resource:22-
biasadd_readvariableop_resource:2
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:22*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџ2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ2
 
_user_specified_nameinputs
Ч	
ѕ
D__inference_dense_8_layer_call_and_return_conditional_losses_1184555

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:2*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџr
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџw
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ2
 
_user_specified_nameinputs
Ч

)__inference_dense_6_layer_call_fn_1184505

inputs
unknown:2
	unknown_0:2
identityЂStatefulPartitionedCallо
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_dense_6_layer_call_and_return_conditional_losses_1183840o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџ2`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
Ч

)__inference_dense_7_layer_call_fn_1184525

inputs
unknown:22
	unknown_0:2
identityЂStatefulPartitionedCallо
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_dense_7_layer_call_and_return_conditional_losses_1183857o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџ2`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ2: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ2
 
_user_specified_nameinputs

ў
I__inference_sequential_2_layer_call_and_return_conditional_losses_1184496

inputs8
&dense_6_matmul_readvariableop_resource:25
'dense_6_biasadd_readvariableop_resource:28
&dense_7_matmul_readvariableop_resource:225
'dense_7_biasadd_readvariableop_resource:28
&dense_8_matmul_readvariableop_resource:25
'dense_8_biasadd_readvariableop_resource:
identityЂdense_6/BiasAdd/ReadVariableOpЂdense_6/MatMul/ReadVariableOpЂdense_7/BiasAdd/ReadVariableOpЂdense_7/MatMul/ReadVariableOpЂdense_8/BiasAdd/ReadVariableOpЂdense_8/MatMul/ReadVariableOp
dense_6/MatMul/ReadVariableOpReadVariableOp&dense_6_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0y
dense_6/MatMulMatMulinputs%dense_6/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dense_6/BiasAdd/ReadVariableOpReadVariableOp'dense_6_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0
dense_6/BiasAddBiasAdddense_6/MatMul:product:0&dense_6/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2`
dense_6/ReluReludense_6/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dense_7/MatMul/ReadVariableOpReadVariableOp&dense_7_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0
dense_7/MatMulMatMuldense_6/Relu:activations:0%dense_7/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dense_7/BiasAdd/ReadVariableOpReadVariableOp'dense_7_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0
dense_7/BiasAddBiasAdddense_7/MatMul:product:0&dense_7/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2`
dense_7/ReluReludense_7/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dense_8/MatMul/ReadVariableOpReadVariableOp&dense_8_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0
dense_8/MatMulMatMuldense_7/Relu:activations:0%dense_8/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
dense_8/BiasAdd/ReadVariableOpReadVariableOp'dense_8_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0
dense_8/BiasAddBiasAdddense_8/MatMul:product:0&dense_8/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџg
IdentityIdentitydense_8/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџ
NoOpNoOp^dense_6/BiasAdd/ReadVariableOp^dense_6/MatMul/ReadVariableOp^dense_7/BiasAdd/ReadVariableOp^dense_7/MatMul/ReadVariableOp^dense_8/BiasAdd/ReadVariableOp^dense_8/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 2@
dense_6/BiasAdd/ReadVariableOpdense_6/BiasAdd/ReadVariableOp2>
dense_6/MatMul/ReadVariableOpdense_6/MatMul/ReadVariableOp2@
dense_7/BiasAdd/ReadVariableOpdense_7/BiasAdd/ReadVariableOp2>
dense_7/MatMul/ReadVariableOpdense_7/MatMul/ReadVariableOp2@
dense_8/BiasAdd/ReadVariableOpdense_8/BiasAdd/ReadVariableOp2>
dense_8/MatMul/ReadVariableOpdense_8/MatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
љ

.__inference_sequential_2_layer_call_fn_1184431

inputs
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityЂStatefulPartitionedCall
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *R
fMRK
I__inference_sequential_2_layer_call_and_return_conditional_losses_1183880o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџ`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
Ч	
ѕ
D__inference_dense_8_layer_call_and_return_conditional_losses_1183873

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:2*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџr
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџw
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ2
 
_user_specified_nameinputs


I__inference_sequential_2_layer_call_and_return_conditional_losses_1184014
dense_6_input!
dense_6_1183998:2
dense_6_1184000:2!
dense_7_1184003:22
dense_7_1184005:2!
dense_8_1184008:2
dense_8_1184010:
identityЂdense_6/StatefulPartitionedCallЂdense_7/StatefulPartitionedCallЂdense_8/StatefulPartitionedCallћ
dense_6/StatefulPartitionedCallStatefulPartitionedCalldense_6_inputdense_6_1183998dense_6_1184000*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_dense_6_layer_call_and_return_conditional_losses_1183840
dense_7/StatefulPartitionedCallStatefulPartitionedCall(dense_6/StatefulPartitionedCall:output:0dense_7_1184003dense_7_1184005*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_dense_7_layer_call_and_return_conditional_losses_1183857
dense_8/StatefulPartitionedCallStatefulPartitionedCall(dense_7/StatefulPartitionedCall:output:0dense_8_1184008dense_8_1184010*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_dense_8_layer_call_and_return_conditional_losses_1183873w
IdentityIdentity(dense_8/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџЌ
NoOpNoOp ^dense_6/StatefulPartitionedCall ^dense_7/StatefulPartitionedCall ^dense_8/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 2B
dense_6/StatefulPartitionedCalldense_6/StatefulPartitionedCall2B
dense_7/StatefulPartitionedCalldense_7/StatefulPartitionedCall2B
dense_8/StatefulPartitionedCalldense_8/StatefulPartitionedCall:V R
'
_output_shapes
:џџџџџџџџџ
'
_user_specified_namedense_6_input
Й
P
$__inference__update_step_xla_1184389
gradient
variable:2*
_XlaMustCompile(*(
_construction_contextkEagerRuntime*
_input_shapes
:2: *
	_noinline(:H D

_output_shapes

:2
"
_user_specified_name
gradient:($
"
_user_specified_name
variable
д
џ
%__inference_signature_wrapper_1184256
input_3
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityЂStatefulPartitionedCallѕ
StatefulPartitionedCallStatefulPartitionedCallinput_3unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *+
f&R$
"__inference__wrapped_model_1183822s
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*+
_output_shapes
:џџџџџџџџџ`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:џџџџџџџџџ
!
_user_specified_name	input_3
J
В
D__inference_model_2_layer_call_and_return_conditional_losses_1184337

inputsE
3sequential_2_dense_6_matmul_readvariableop_resource:2B
4sequential_2_dense_6_biasadd_readvariableop_resource:2E
3sequential_2_dense_7_matmul_readvariableop_resource:22B
4sequential_2_dense_7_biasadd_readvariableop_resource:2E
3sequential_2_dense_8_matmul_readvariableop_resource:2B
4sequential_2_dense_8_biasadd_readvariableop_resource:
identityЂ+sequential_2/dense_6/BiasAdd/ReadVariableOpЂ-sequential_2/dense_6/BiasAdd_1/ReadVariableOpЂ*sequential_2/dense_6/MatMul/ReadVariableOpЂ,sequential_2/dense_6/MatMul_1/ReadVariableOpЂ+sequential_2/dense_7/BiasAdd/ReadVariableOpЂ-sequential_2/dense_7/BiasAdd_1/ReadVariableOpЂ*sequential_2/dense_7/MatMul/ReadVariableOpЂ,sequential_2/dense_7/MatMul_1/ReadVariableOpЂ+sequential_2/dense_8/BiasAdd/ReadVariableOpЂ-sequential_2/dense_8/BiasAdd_1/ReadVariableOpЂ*sequential_2/dense_8/MatMul/ReadVariableOpЂ,sequential_2/dense_8/MatMul_1/ReadVariableOp
.tf.__operators__.getitem_5/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       
0tf.__operators__.getitem_5/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        
0tf.__operators__.getitem_5/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      а
(tf.__operators__.getitem_5/strided_sliceStridedSliceinputs7tf.__operators__.getitem_5/strided_slice/stack:output:09tf.__operators__.getitem_5/strided_slice/stack_1:output:09tf.__operators__.getitem_5/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*

begin_mask*
end_mask
.tf.__operators__.getitem_4/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        
0tf.__operators__.getitem_4/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       
0tf.__operators__.getitem_4/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      а
(tf.__operators__.getitem_4/strided_sliceStridedSliceinputs7tf.__operators__.getitem_4/strided_slice/stack:output:09tf.__operators__.getitem_4/strided_slice/stack_1:output:09tf.__operators__.getitem_4/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*

begin_mask*
end_mask
*sequential_2/dense_6/MatMul/ReadVariableOpReadVariableOp3sequential_2_dense_6_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0О
sequential_2/dense_6/MatMulMatMul1tf.__operators__.getitem_4/strided_slice:output:02sequential_2/dense_6/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
+sequential_2/dense_6/BiasAdd/ReadVariableOpReadVariableOp4sequential_2_dense_6_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0Е
sequential_2/dense_6/BiasAddBiasAdd%sequential_2/dense_6/MatMul:product:03sequential_2/dense_6/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2z
sequential_2/dense_6/ReluRelu%sequential_2/dense_6/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
*sequential_2/dense_7/MatMul/ReadVariableOpReadVariableOp3sequential_2_dense_7_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Д
sequential_2/dense_7/MatMulMatMul'sequential_2/dense_6/Relu:activations:02sequential_2/dense_7/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
+sequential_2/dense_7/BiasAdd/ReadVariableOpReadVariableOp4sequential_2_dense_7_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0Е
sequential_2/dense_7/BiasAddBiasAdd%sequential_2/dense_7/MatMul:product:03sequential_2/dense_7/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2z
sequential_2/dense_7/ReluRelu%sequential_2/dense_7/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
*sequential_2/dense_8/MatMul/ReadVariableOpReadVariableOp3sequential_2_dense_8_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0Д
sequential_2/dense_8/MatMulMatMul'sequential_2/dense_7/Relu:activations:02sequential_2/dense_8/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
+sequential_2/dense_8/BiasAdd/ReadVariableOpReadVariableOp4sequential_2_dense_8_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0Е
sequential_2/dense_8/BiasAddBiasAdd%sequential_2/dense_8/MatMul:product:03sequential_2/dense_8/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ 
,sequential_2/dense_6/MatMul_1/ReadVariableOpReadVariableOp3sequential_2_dense_6_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0Т
sequential_2/dense_6/MatMul_1MatMul1tf.__operators__.getitem_5/strided_slice:output:04sequential_2/dense_6/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
-sequential_2/dense_6/BiasAdd_1/ReadVariableOpReadVariableOp4sequential_2_dense_6_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0Л
sequential_2/dense_6/BiasAdd_1BiasAdd'sequential_2/dense_6/MatMul_1:product:05sequential_2/dense_6/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2~
sequential_2/dense_6/Relu_1Relu'sequential_2/dense_6/BiasAdd_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
,sequential_2/dense_7/MatMul_1/ReadVariableOpReadVariableOp3sequential_2_dense_7_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0К
sequential_2/dense_7/MatMul_1MatMul)sequential_2/dense_6/Relu_1:activations:04sequential_2/dense_7/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
-sequential_2/dense_7/BiasAdd_1/ReadVariableOpReadVariableOp4sequential_2_dense_7_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0Л
sequential_2/dense_7/BiasAdd_1BiasAdd'sequential_2/dense_7/MatMul_1:product:05sequential_2/dense_7/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2~
sequential_2/dense_7/Relu_1Relu'sequential_2/dense_7/BiasAdd_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
,sequential_2/dense_8/MatMul_1/ReadVariableOpReadVariableOp3sequential_2_dense_8_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0К
sequential_2/dense_8/MatMul_1MatMul)sequential_2/dense_7/Relu_1:activations:04sequential_2/dense_8/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
-sequential_2/dense_8/BiasAdd_1/ReadVariableOpReadVariableOp4sequential_2_dense_8_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0Л
sequential_2/dense_8/BiasAdd_1BiasAdd'sequential_2/dense_8/MatMul_1:product:05sequential_2/dense_8/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџГ
tf.stack_2/stackPack%sequential_2/dense_8/BiasAdd:output:0'sequential_2/dense_8/BiasAdd_1:output:0*
N*
T0*+
_output_shapes
:џџџџџџџџџ*

axisl
IdentityIdentitytf.stack_2/stack:output:0^NoOp*
T0*+
_output_shapes
:џџџџџџџџџє
NoOpNoOp,^sequential_2/dense_6/BiasAdd/ReadVariableOp.^sequential_2/dense_6/BiasAdd_1/ReadVariableOp+^sequential_2/dense_6/MatMul/ReadVariableOp-^sequential_2/dense_6/MatMul_1/ReadVariableOp,^sequential_2/dense_7/BiasAdd/ReadVariableOp.^sequential_2/dense_7/BiasAdd_1/ReadVariableOp+^sequential_2/dense_7/MatMul/ReadVariableOp-^sequential_2/dense_7/MatMul_1/ReadVariableOp,^sequential_2/dense_8/BiasAdd/ReadVariableOp.^sequential_2/dense_8/BiasAdd_1/ReadVariableOp+^sequential_2/dense_8/MatMul/ReadVariableOp-^sequential_2/dense_8/MatMul_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 2Z
+sequential_2/dense_6/BiasAdd/ReadVariableOp+sequential_2/dense_6/BiasAdd/ReadVariableOp2^
-sequential_2/dense_6/BiasAdd_1/ReadVariableOp-sequential_2/dense_6/BiasAdd_1/ReadVariableOp2X
*sequential_2/dense_6/MatMul/ReadVariableOp*sequential_2/dense_6/MatMul/ReadVariableOp2\
,sequential_2/dense_6/MatMul_1/ReadVariableOp,sequential_2/dense_6/MatMul_1/ReadVariableOp2Z
+sequential_2/dense_7/BiasAdd/ReadVariableOp+sequential_2/dense_7/BiasAdd/ReadVariableOp2^
-sequential_2/dense_7/BiasAdd_1/ReadVariableOp-sequential_2/dense_7/BiasAdd_1/ReadVariableOp2X
*sequential_2/dense_7/MatMul/ReadVariableOp*sequential_2/dense_7/MatMul/ReadVariableOp2\
,sequential_2/dense_7/MatMul_1/ReadVariableOp,sequential_2/dense_7/MatMul_1/ReadVariableOp2Z
+sequential_2/dense_8/BiasAdd/ReadVariableOp+sequential_2/dense_8/BiasAdd/ReadVariableOp2^
-sequential_2/dense_8/BiasAdd_1/ReadVariableOp-sequential_2/dense_8/BiasAdd_1/ReadVariableOp2X
*sequential_2/dense_8/MatMul/ReadVariableOp*sequential_2/dense_8/MatMul/ReadVariableOp2\
,sequential_2/dense_8/MatMul_1/ReadVariableOp,sequential_2/dense_8/MatMul_1/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs

ў
I__inference_sequential_2_layer_call_and_return_conditional_losses_1184472

inputs8
&dense_6_matmul_readvariableop_resource:25
'dense_6_biasadd_readvariableop_resource:28
&dense_7_matmul_readvariableop_resource:225
'dense_7_biasadd_readvariableop_resource:28
&dense_8_matmul_readvariableop_resource:25
'dense_8_biasadd_readvariableop_resource:
identityЂdense_6/BiasAdd/ReadVariableOpЂdense_6/MatMul/ReadVariableOpЂdense_7/BiasAdd/ReadVariableOpЂdense_7/MatMul/ReadVariableOpЂdense_8/BiasAdd/ReadVariableOpЂdense_8/MatMul/ReadVariableOp
dense_6/MatMul/ReadVariableOpReadVariableOp&dense_6_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0y
dense_6/MatMulMatMulinputs%dense_6/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dense_6/BiasAdd/ReadVariableOpReadVariableOp'dense_6_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0
dense_6/BiasAddBiasAdddense_6/MatMul:product:0&dense_6/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2`
dense_6/ReluReludense_6/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dense_7/MatMul/ReadVariableOpReadVariableOp&dense_7_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0
dense_7/MatMulMatMuldense_6/Relu:activations:0%dense_7/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dense_7/BiasAdd/ReadVariableOpReadVariableOp'dense_7_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0
dense_7/BiasAddBiasAdddense_7/MatMul:product:0&dense_7/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2`
dense_7/ReluReludense_7/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dense_8/MatMul/ReadVariableOpReadVariableOp&dense_8_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0
dense_8/MatMulMatMuldense_7/Relu:activations:0%dense_8/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
dense_8/BiasAdd/ReadVariableOpReadVariableOp'dense_8_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0
dense_8/BiasAddBiasAdddense_8/MatMul:product:0&dense_8/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџg
IdentityIdentitydense_8/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџ
NoOpNoOp^dense_6/BiasAdd/ReadVariableOp^dense_6/MatMul/ReadVariableOp^dense_7/BiasAdd/ReadVariableOp^dense_7/MatMul/ReadVariableOp^dense_8/BiasAdd/ReadVariableOp^dense_8/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 2@
dense_6/BiasAdd/ReadVariableOpdense_6/BiasAdd/ReadVariableOp2>
dense_6/MatMul/ReadVariableOpdense_6/MatMul/ReadVariableOp2@
dense_7/BiasAdd/ReadVariableOpdense_7/BiasAdd/ReadVariableOp2>
dense_7/MatMul/ReadVariableOpdense_7/MatMul/ReadVariableOp2@
dense_8/BiasAdd/ReadVariableOpdense_8/BiasAdd/ReadVariableOp2>
dense_8/MatMul/ReadVariableOpdense_8/MatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
П

D__inference_model_2_layer_call_and_return_conditional_losses_1184202
input_3&
sequential_2_1184180:2"
sequential_2_1184182:2&
sequential_2_1184184:22"
sequential_2_1184186:2&
sequential_2_1184188:2"
sequential_2_1184190:
identityЂ$sequential_2/StatefulPartitionedCallЂ&sequential_2/StatefulPartitionedCall_1
.tf.__operators__.getitem_5/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       
0tf.__operators__.getitem_5/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        
0tf.__operators__.getitem_5/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      б
(tf.__operators__.getitem_5/strided_sliceStridedSliceinput_37tf.__operators__.getitem_5/strided_slice/stack:output:09tf.__operators__.getitem_5/strided_slice/stack_1:output:09tf.__operators__.getitem_5/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*

begin_mask*
end_mask
.tf.__operators__.getitem_4/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        
0tf.__operators__.getitem_4/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       
0tf.__operators__.getitem_4/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      б
(tf.__operators__.getitem_4/strided_sliceStridedSliceinput_37tf.__operators__.getitem_4/strided_slice/stack:output:09tf.__operators__.getitem_4/strided_slice/stack_1:output:09tf.__operators__.getitem_4/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*

begin_mask*
end_mask
$sequential_2/StatefulPartitionedCallStatefulPartitionedCall1tf.__operators__.getitem_4/strided_slice:output:0sequential_2_1184180sequential_2_1184182sequential_2_1184184sequential_2_1184186sequential_2_1184188sequential_2_1184190*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *R
fMRK
I__inference_sequential_2_layer_call_and_return_conditional_losses_1183880
&sequential_2/StatefulPartitionedCall_1StatefulPartitionedCall1tf.__operators__.getitem_5/strided_slice:output:0sequential_2_1184180sequential_2_1184182sequential_2_1184184sequential_2_1184186sequential_2_1184188sequential_2_1184190*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *R
fMRK
I__inference_sequential_2_layer_call_and_return_conditional_losses_1183880У
tf.stack_2/stackPack-sequential_2/StatefulPartitionedCall:output:0/sequential_2/StatefulPartitionedCall_1:output:0*
N*
T0*+
_output_shapes
:џџџџџџџџџ*

axisl
IdentityIdentitytf.stack_2/stack:output:0^NoOp*
T0*+
_output_shapes
:џџџџџџџџџ
NoOpNoOp%^sequential_2/StatefulPartitionedCall'^sequential_2/StatefulPartitionedCall_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 2L
$sequential_2/StatefulPartitionedCall$sequential_2/StatefulPartitionedCall2P
&sequential_2/StatefulPartitionedCall_1&sequential_2/StatefulPartitionedCall_1:P L
'
_output_shapes
:џџџџџџџџџ
!
_user_specified_name	input_3


ѕ
D__inference_dense_6_layer_call_and_return_conditional_losses_1184516

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:2
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:2*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџ2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ї

)__inference_model_2_layer_call_fn_1184290

inputs
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityЂStatefulPartitionedCall
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_model_2_layer_call_and_return_conditional_losses_1184137s
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*+
_output_shapes
:џџџџџџџџџ`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
П

D__inference_model_2_layer_call_and_return_conditional_losses_1184235
input_3&
sequential_2_1184213:2"
sequential_2_1184215:2&
sequential_2_1184217:22"
sequential_2_1184219:2&
sequential_2_1184221:2"
sequential_2_1184223:
identityЂ$sequential_2/StatefulPartitionedCallЂ&sequential_2/StatefulPartitionedCall_1
.tf.__operators__.getitem_5/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       
0tf.__operators__.getitem_5/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        
0tf.__operators__.getitem_5/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      б
(tf.__operators__.getitem_5/strided_sliceStridedSliceinput_37tf.__operators__.getitem_5/strided_slice/stack:output:09tf.__operators__.getitem_5/strided_slice/stack_1:output:09tf.__operators__.getitem_5/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*

begin_mask*
end_mask
.tf.__operators__.getitem_4/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        
0tf.__operators__.getitem_4/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       
0tf.__operators__.getitem_4/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      б
(tf.__operators__.getitem_4/strided_sliceStridedSliceinput_37tf.__operators__.getitem_4/strided_slice/stack:output:09tf.__operators__.getitem_4/strided_slice/stack_1:output:09tf.__operators__.getitem_4/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*

begin_mask*
end_mask
$sequential_2/StatefulPartitionedCallStatefulPartitionedCall1tf.__operators__.getitem_4/strided_slice:output:0sequential_2_1184213sequential_2_1184215sequential_2_1184217sequential_2_1184219sequential_2_1184221sequential_2_1184223*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *R
fMRK
I__inference_sequential_2_layer_call_and_return_conditional_losses_1183963
&sequential_2/StatefulPartitionedCall_1StatefulPartitionedCall1tf.__operators__.getitem_5/strided_slice:output:0sequential_2_1184213sequential_2_1184215sequential_2_1184217sequential_2_1184219sequential_2_1184221sequential_2_1184223*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *R
fMRK
I__inference_sequential_2_layer_call_and_return_conditional_losses_1183963У
tf.stack_2/stackPack-sequential_2/StatefulPartitionedCall:output:0/sequential_2/StatefulPartitionedCall_1:output:0*
N*
T0*+
_output_shapes
:џџџџџџџџџ*

axisl
IdentityIdentitytf.stack_2/stack:output:0^NoOp*
T0*+
_output_shapes
:џџџџџџџџџ
NoOpNoOp%^sequential_2/StatefulPartitionedCall'^sequential_2/StatefulPartitionedCall_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 2L
$sequential_2/StatefulPartitionedCall$sequential_2/StatefulPartitionedCall2P
&sequential_2/StatefulPartitionedCall_1&sequential_2/StatefulPartitionedCall_1:P L
'
_output_shapes
:џџџџџџџџџ
!
_user_specified_name	input_3
Й
P
$__inference__update_step_xla_1184409
gradient
variable:2*
_XlaMustCompile(*(
_construction_contextkEagerRuntime*
_input_shapes
:2: *
	_noinline(:H D

_output_shapes

:2
"
_user_specified_name
gradient:($
"
_user_specified_name
variable
њ

)__inference_model_2_layer_call_fn_1184085
input_3
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityЂStatefulPartitionedCall
StatefulPartitionedCallStatefulPartitionedCallinput_3unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_model_2_layer_call_and_return_conditional_losses_1184070s
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*+
_output_shapes
:џџџџџџџџџ`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:џџџџџџџџџ
!
_user_specified_name	input_3
­
L
$__inference__update_step_xla_1184414
gradient
variable:*
_XlaMustCompile(*(
_construction_contextkEagerRuntime*
_input_shapes

:: *
	_noinline(:D @

_output_shapes
:
"
_user_specified_name
gradient:($
"
_user_specified_name
variable
ЮR
Ё	
"__inference__wrapped_model_1183822
input_3M
;model_2_sequential_2_dense_6_matmul_readvariableop_resource:2J
<model_2_sequential_2_dense_6_biasadd_readvariableop_resource:2M
;model_2_sequential_2_dense_7_matmul_readvariableop_resource:22J
<model_2_sequential_2_dense_7_biasadd_readvariableop_resource:2M
;model_2_sequential_2_dense_8_matmul_readvariableop_resource:2J
<model_2_sequential_2_dense_8_biasadd_readvariableop_resource:
identityЂ3model_2/sequential_2/dense_6/BiasAdd/ReadVariableOpЂ5model_2/sequential_2/dense_6/BiasAdd_1/ReadVariableOpЂ2model_2/sequential_2/dense_6/MatMul/ReadVariableOpЂ4model_2/sequential_2/dense_6/MatMul_1/ReadVariableOpЂ3model_2/sequential_2/dense_7/BiasAdd/ReadVariableOpЂ5model_2/sequential_2/dense_7/BiasAdd_1/ReadVariableOpЂ2model_2/sequential_2/dense_7/MatMul/ReadVariableOpЂ4model_2/sequential_2/dense_7/MatMul_1/ReadVariableOpЂ3model_2/sequential_2/dense_8/BiasAdd/ReadVariableOpЂ5model_2/sequential_2/dense_8/BiasAdd_1/ReadVariableOpЂ2model_2/sequential_2/dense_8/MatMul/ReadVariableOpЂ4model_2/sequential_2/dense_8/MatMul_1/ReadVariableOp
6model_2/tf.__operators__.getitem_5/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       
8model_2/tf.__operators__.getitem_5/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        
8model_2/tf.__operators__.getitem_5/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ё
0model_2/tf.__operators__.getitem_5/strided_sliceStridedSliceinput_3?model_2/tf.__operators__.getitem_5/strided_slice/stack:output:0Amodel_2/tf.__operators__.getitem_5/strided_slice/stack_1:output:0Amodel_2/tf.__operators__.getitem_5/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*

begin_mask*
end_mask
6model_2/tf.__operators__.getitem_4/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        
8model_2/tf.__operators__.getitem_4/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       
8model_2/tf.__operators__.getitem_4/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ё
0model_2/tf.__operators__.getitem_4/strided_sliceStridedSliceinput_3?model_2/tf.__operators__.getitem_4/strided_slice/stack:output:0Amodel_2/tf.__operators__.getitem_4/strided_slice/stack_1:output:0Amodel_2/tf.__operators__.getitem_4/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*

begin_mask*
end_maskЎ
2model_2/sequential_2/dense_6/MatMul/ReadVariableOpReadVariableOp;model_2_sequential_2_dense_6_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0ж
#model_2/sequential_2/dense_6/MatMulMatMul9model_2/tf.__operators__.getitem_4/strided_slice:output:0:model_2/sequential_2/dense_6/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2Ќ
3model_2/sequential_2/dense_6/BiasAdd/ReadVariableOpReadVariableOp<model_2_sequential_2_dense_6_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0Э
$model_2/sequential_2/dense_6/BiasAddBiasAdd-model_2/sequential_2/dense_6/MatMul:product:0;model_2/sequential_2/dense_6/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
!model_2/sequential_2/dense_6/ReluRelu-model_2/sequential_2/dense_6/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2Ў
2model_2/sequential_2/dense_7/MatMul/ReadVariableOpReadVariableOp;model_2_sequential_2_dense_7_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Ь
#model_2/sequential_2/dense_7/MatMulMatMul/model_2/sequential_2/dense_6/Relu:activations:0:model_2/sequential_2/dense_7/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2Ќ
3model_2/sequential_2/dense_7/BiasAdd/ReadVariableOpReadVariableOp<model_2_sequential_2_dense_7_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0Э
$model_2/sequential_2/dense_7/BiasAddBiasAdd-model_2/sequential_2/dense_7/MatMul:product:0;model_2/sequential_2/dense_7/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
!model_2/sequential_2/dense_7/ReluRelu-model_2/sequential_2/dense_7/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2Ў
2model_2/sequential_2/dense_8/MatMul/ReadVariableOpReadVariableOp;model_2_sequential_2_dense_8_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0Ь
#model_2/sequential_2/dense_8/MatMulMatMul/model_2/sequential_2/dense_7/Relu:activations:0:model_2/sequential_2/dense_8/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџЌ
3model_2/sequential_2/dense_8/BiasAdd/ReadVariableOpReadVariableOp<model_2_sequential_2_dense_8_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0Э
$model_2/sequential_2/dense_8/BiasAddBiasAdd-model_2/sequential_2/dense_8/MatMul:product:0;model_2/sequential_2/dense_8/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџА
4model_2/sequential_2/dense_6/MatMul_1/ReadVariableOpReadVariableOp;model_2_sequential_2_dense_6_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0к
%model_2/sequential_2/dense_6/MatMul_1MatMul9model_2/tf.__operators__.getitem_5/strided_slice:output:0<model_2/sequential_2/dense_6/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2Ў
5model_2/sequential_2/dense_6/BiasAdd_1/ReadVariableOpReadVariableOp<model_2_sequential_2_dense_6_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0г
&model_2/sequential_2/dense_6/BiasAdd_1BiasAdd/model_2/sequential_2/dense_6/MatMul_1:product:0=model_2/sequential_2/dense_6/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
#model_2/sequential_2/dense_6/Relu_1Relu/model_2/sequential_2/dense_6/BiasAdd_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2А
4model_2/sequential_2/dense_7/MatMul_1/ReadVariableOpReadVariableOp;model_2_sequential_2_dense_7_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0в
%model_2/sequential_2/dense_7/MatMul_1MatMul1model_2/sequential_2/dense_6/Relu_1:activations:0<model_2/sequential_2/dense_7/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2Ў
5model_2/sequential_2/dense_7/BiasAdd_1/ReadVariableOpReadVariableOp<model_2_sequential_2_dense_7_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0г
&model_2/sequential_2/dense_7/BiasAdd_1BiasAdd/model_2/sequential_2/dense_7/MatMul_1:product:0=model_2/sequential_2/dense_7/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
#model_2/sequential_2/dense_7/Relu_1Relu/model_2/sequential_2/dense_7/BiasAdd_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2А
4model_2/sequential_2/dense_8/MatMul_1/ReadVariableOpReadVariableOp;model_2_sequential_2_dense_8_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0в
%model_2/sequential_2/dense_8/MatMul_1MatMul1model_2/sequential_2/dense_7/Relu_1:activations:0<model_2/sequential_2/dense_8/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџЎ
5model_2/sequential_2/dense_8/BiasAdd_1/ReadVariableOpReadVariableOp<model_2_sequential_2_dense_8_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0г
&model_2/sequential_2/dense_8/BiasAdd_1BiasAdd/model_2/sequential_2/dense_8/MatMul_1:product:0=model_2/sequential_2/dense_8/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџЫ
model_2/tf.stack_2/stackPack-model_2/sequential_2/dense_8/BiasAdd:output:0/model_2/sequential_2/dense_8/BiasAdd_1:output:0*
N*
T0*+
_output_shapes
:џџџџџџџџџ*

axist
IdentityIdentity!model_2/tf.stack_2/stack:output:0^NoOp*
T0*+
_output_shapes
:џџџџџџџџџд
NoOpNoOp4^model_2/sequential_2/dense_6/BiasAdd/ReadVariableOp6^model_2/sequential_2/dense_6/BiasAdd_1/ReadVariableOp3^model_2/sequential_2/dense_6/MatMul/ReadVariableOp5^model_2/sequential_2/dense_6/MatMul_1/ReadVariableOp4^model_2/sequential_2/dense_7/BiasAdd/ReadVariableOp6^model_2/sequential_2/dense_7/BiasAdd_1/ReadVariableOp3^model_2/sequential_2/dense_7/MatMul/ReadVariableOp5^model_2/sequential_2/dense_7/MatMul_1/ReadVariableOp4^model_2/sequential_2/dense_8/BiasAdd/ReadVariableOp6^model_2/sequential_2/dense_8/BiasAdd_1/ReadVariableOp3^model_2/sequential_2/dense_8/MatMul/ReadVariableOp5^model_2/sequential_2/dense_8/MatMul_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 2j
3model_2/sequential_2/dense_6/BiasAdd/ReadVariableOp3model_2/sequential_2/dense_6/BiasAdd/ReadVariableOp2n
5model_2/sequential_2/dense_6/BiasAdd_1/ReadVariableOp5model_2/sequential_2/dense_6/BiasAdd_1/ReadVariableOp2h
2model_2/sequential_2/dense_6/MatMul/ReadVariableOp2model_2/sequential_2/dense_6/MatMul/ReadVariableOp2l
4model_2/sequential_2/dense_6/MatMul_1/ReadVariableOp4model_2/sequential_2/dense_6/MatMul_1/ReadVariableOp2j
3model_2/sequential_2/dense_7/BiasAdd/ReadVariableOp3model_2/sequential_2/dense_7/BiasAdd/ReadVariableOp2n
5model_2/sequential_2/dense_7/BiasAdd_1/ReadVariableOp5model_2/sequential_2/dense_7/BiasAdd_1/ReadVariableOp2h
2model_2/sequential_2/dense_7/MatMul/ReadVariableOp2model_2/sequential_2/dense_7/MatMul/ReadVariableOp2l
4model_2/sequential_2/dense_7/MatMul_1/ReadVariableOp4model_2/sequential_2/dense_7/MatMul_1/ReadVariableOp2j
3model_2/sequential_2/dense_8/BiasAdd/ReadVariableOp3model_2/sequential_2/dense_8/BiasAdd/ReadVariableOp2n
5model_2/sequential_2/dense_8/BiasAdd_1/ReadVariableOp5model_2/sequential_2/dense_8/BiasAdd_1/ReadVariableOp2h
2model_2/sequential_2/dense_8/MatMul/ReadVariableOp2model_2/sequential_2/dense_8/MatMul/ReadVariableOp2l
4model_2/sequential_2/dense_8/MatMul_1/ReadVariableOp4model_2/sequential_2/dense_8/MatMul_1/ReadVariableOp:P L
'
_output_shapes
:џџџџџџџџџ
!
_user_specified_name	input_3
юe

#__inference__traced_restore_1184732
file_prefix1
assignvariableop_dense_6_kernel:2-
assignvariableop_1_dense_6_bias:23
!assignvariableop_2_dense_7_kernel:22-
assignvariableop_3_dense_7_bias:23
!assignvariableop_4_dense_8_kernel:2-
assignvariableop_5_dense_8_bias:&
assignvariableop_6_iteration:	 *
 assignvariableop_7_learning_rate: :
(assignvariableop_8_adam_m_dense_6_kernel:2:
(assignvariableop_9_adam_v_dense_6_kernel:25
'assignvariableop_10_adam_m_dense_6_bias:25
'assignvariableop_11_adam_v_dense_6_bias:2;
)assignvariableop_12_adam_m_dense_7_kernel:22;
)assignvariableop_13_adam_v_dense_7_kernel:225
'assignvariableop_14_adam_m_dense_7_bias:25
'assignvariableop_15_adam_v_dense_7_bias:2;
)assignvariableop_16_adam_m_dense_8_kernel:2;
)assignvariableop_17_adam_v_dense_8_kernel:25
'assignvariableop_18_adam_m_dense_8_bias:5
'assignvariableop_19_adam_v_dense_8_bias:%
assignvariableop_20_total_1: %
assignvariableop_21_count_1: #
assignvariableop_22_total: #
assignvariableop_23_count: 
identity_25ЂAssignVariableOpЂAssignVariableOp_1ЂAssignVariableOp_10ЂAssignVariableOp_11ЂAssignVariableOp_12ЂAssignVariableOp_13ЂAssignVariableOp_14ЂAssignVariableOp_15ЂAssignVariableOp_16ЂAssignVariableOp_17ЂAssignVariableOp_18ЂAssignVariableOp_19ЂAssignVariableOp_2ЂAssignVariableOp_20ЂAssignVariableOp_21ЂAssignVariableOp_22ЂAssignVariableOp_23ЂAssignVariableOp_3ЂAssignVariableOp_4ЂAssignVariableOp_5ЂAssignVariableOp_6ЂAssignVariableOp_7ЂAssignVariableOp_8ЂAssignVariableOp_9Ѓ

RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*Щ	
valueП	BМ	B&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUEB3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPHЂ
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*E
value<B:B B B B B B B B B B B B B B B B B B B B B B B B B 
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*x
_output_shapesf
d:::::::::::::::::::::::::*'
dtypes
2	[
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:В
AssignVariableOpAssignVariableOpassignvariableop_dense_6_kernelIdentity:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:Ж
AssignVariableOp_1AssignVariableOpassignvariableop_1_dense_6_biasIdentity_1:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:И
AssignVariableOp_2AssignVariableOp!assignvariableop_2_dense_7_kernelIdentity_2:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:Ж
AssignVariableOp_3AssignVariableOpassignvariableop_3_dense_7_biasIdentity_3:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:И
AssignVariableOp_4AssignVariableOp!assignvariableop_4_dense_8_kernelIdentity_4:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:Ж
AssignVariableOp_5AssignVariableOpassignvariableop_5_dense_8_biasIdentity_5:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0	*
_output_shapes
:Г
AssignVariableOp_6AssignVariableOpassignvariableop_6_iterationIdentity_6:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0	]

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:З
AssignVariableOp_7AssignVariableOp assignvariableop_7_learning_rateIdentity_7:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:П
AssignVariableOp_8AssignVariableOp(assignvariableop_8_adam_m_dense_6_kernelIdentity_8:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:П
AssignVariableOp_9AssignVariableOp(assignvariableop_9_adam_v_dense_6_kernelIdentity_9:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:Р
AssignVariableOp_10AssignVariableOp'assignvariableop_10_adam_m_dense_6_biasIdentity_10:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:Р
AssignVariableOp_11AssignVariableOp'assignvariableop_11_adam_v_dense_6_biasIdentity_11:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:Т
AssignVariableOp_12AssignVariableOp)assignvariableop_12_adam_m_dense_7_kernelIdentity_12:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:Т
AssignVariableOp_13AssignVariableOp)assignvariableop_13_adam_v_dense_7_kernelIdentity_13:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:Р
AssignVariableOp_14AssignVariableOp'assignvariableop_14_adam_m_dense_7_biasIdentity_14:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:Р
AssignVariableOp_15AssignVariableOp'assignvariableop_15_adam_v_dense_7_biasIdentity_15:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:Т
AssignVariableOp_16AssignVariableOp)assignvariableop_16_adam_m_dense_8_kernelIdentity_16:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:Т
AssignVariableOp_17AssignVariableOp)assignvariableop_17_adam_v_dense_8_kernelIdentity_17:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:Р
AssignVariableOp_18AssignVariableOp'assignvariableop_18_adam_m_dense_8_biasIdentity_18:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:Р
AssignVariableOp_19AssignVariableOp'assignvariableop_19_adam_v_dense_8_biasIdentity_19:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:Д
AssignVariableOp_20AssignVariableOpassignvariableop_20_total_1Identity_20:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:Д
AssignVariableOp_21AssignVariableOpassignvariableop_21_count_1Identity_21:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:В
AssignVariableOp_22AssignVariableOpassignvariableop_22_totalIdentity_22:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:В
AssignVariableOp_23AssignVariableOpassignvariableop_23_countIdentity_23:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0Y
NoOpNoOp"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 п
Identity_24Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: W
Identity_25IdentityIdentity_24:output:0^NoOp_1*
T0*
_output_shapes
: Ь
NoOp_1NoOp^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9*"
_acd_function_control_output(*
_output_shapes
 "#
identity_25Identity_25:output:0*E
_input_shapes4
2: : : : : : : : : : : : : : : : : : : : : : : : : 2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142*
AssignVariableOp_15AssignVariableOp_152*
AssignVariableOp_16AssignVariableOp_162*
AssignVariableOp_17AssignVariableOp_172*
AssignVariableOp_18AssignVariableOp_182*
AssignVariableOp_19AssignVariableOp_192(
AssignVariableOp_2AssignVariableOp_22*
AssignVariableOp_20AssignVariableOp_202*
AssignVariableOp_21AssignVariableOp_212*
AssignVariableOp_22AssignVariableOp_222*
AssignVariableOp_23AssignVariableOp_232(
AssignVariableOp_3AssignVariableOp_32(
AssignVariableOp_4AssignVariableOp_42(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_9:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
	

.__inference_sequential_2_layer_call_fn_1183895
dense_6_input
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityЂStatefulPartitionedCall
StatefulPartitionedCallStatefulPartitionedCalldense_6_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *R
fMRK
I__inference_sequential_2_layer_call_and_return_conditional_losses_1183880o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџ`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:V R
'
_output_shapes
:џџџџџџџџџ
'
_user_specified_namedense_6_input
Л

D__inference_model_2_layer_call_and_return_conditional_losses_1184070

inputs&
sequential_2_1184048:2"
sequential_2_1184050:2&
sequential_2_1184052:22"
sequential_2_1184054:2&
sequential_2_1184056:2"
sequential_2_1184058:
identityЂ$sequential_2/StatefulPartitionedCallЂ&sequential_2/StatefulPartitionedCall_1
.tf.__operators__.getitem_5/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       
0tf.__operators__.getitem_5/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        
0tf.__operators__.getitem_5/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      а
(tf.__operators__.getitem_5/strided_sliceStridedSliceinputs7tf.__operators__.getitem_5/strided_slice/stack:output:09tf.__operators__.getitem_5/strided_slice/stack_1:output:09tf.__operators__.getitem_5/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*

begin_mask*
end_mask
.tf.__operators__.getitem_4/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        
0tf.__operators__.getitem_4/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       
0tf.__operators__.getitem_4/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      а
(tf.__operators__.getitem_4/strided_sliceStridedSliceinputs7tf.__operators__.getitem_4/strided_slice/stack:output:09tf.__operators__.getitem_4/strided_slice/stack_1:output:09tf.__operators__.getitem_4/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*

begin_mask*
end_mask
$sequential_2/StatefulPartitionedCallStatefulPartitionedCall1tf.__operators__.getitem_4/strided_slice:output:0sequential_2_1184048sequential_2_1184050sequential_2_1184052sequential_2_1184054sequential_2_1184056sequential_2_1184058*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *R
fMRK
I__inference_sequential_2_layer_call_and_return_conditional_losses_1183880
&sequential_2/StatefulPartitionedCall_1StatefulPartitionedCall1tf.__operators__.getitem_5/strided_slice:output:0sequential_2_1184048sequential_2_1184050sequential_2_1184052sequential_2_1184054sequential_2_1184056sequential_2_1184058*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *R
fMRK
I__inference_sequential_2_layer_call_and_return_conditional_losses_1183880У
tf.stack_2/stackPack-sequential_2/StatefulPartitionedCall:output:0/sequential_2/StatefulPartitionedCall_1:output:0*
N*
T0*+
_output_shapes
:џџџџџџџџџ*

axisl
IdentityIdentitytf.stack_2/stack:output:0^NoOp*
T0*+
_output_shapes
:џџџџџџџџџ
NoOpNoOp%^sequential_2/StatefulPartitionedCall'^sequential_2/StatefulPartitionedCall_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 2L
$sequential_2/StatefulPartitionedCall$sequential_2/StatefulPartitionedCall2P
&sequential_2/StatefulPartitionedCall_1&sequential_2/StatefulPartitionedCall_1:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
ї

)__inference_model_2_layer_call_fn_1184273

inputs
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityЂStatefulPartitionedCall
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_model_2_layer_call_and_return_conditional_losses_1184070s
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*+
_output_shapes
:џџџџџџџџџ`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs


ѕ
D__inference_dense_6_layer_call_and_return_conditional_losses_1183840

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:2
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:2*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџ2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
J
В
D__inference_model_2_layer_call_and_return_conditional_losses_1184384

inputsE
3sequential_2_dense_6_matmul_readvariableop_resource:2B
4sequential_2_dense_6_biasadd_readvariableop_resource:2E
3sequential_2_dense_7_matmul_readvariableop_resource:22B
4sequential_2_dense_7_biasadd_readvariableop_resource:2E
3sequential_2_dense_8_matmul_readvariableop_resource:2B
4sequential_2_dense_8_biasadd_readvariableop_resource:
identityЂ+sequential_2/dense_6/BiasAdd/ReadVariableOpЂ-sequential_2/dense_6/BiasAdd_1/ReadVariableOpЂ*sequential_2/dense_6/MatMul/ReadVariableOpЂ,sequential_2/dense_6/MatMul_1/ReadVariableOpЂ+sequential_2/dense_7/BiasAdd/ReadVariableOpЂ-sequential_2/dense_7/BiasAdd_1/ReadVariableOpЂ*sequential_2/dense_7/MatMul/ReadVariableOpЂ,sequential_2/dense_7/MatMul_1/ReadVariableOpЂ+sequential_2/dense_8/BiasAdd/ReadVariableOpЂ-sequential_2/dense_8/BiasAdd_1/ReadVariableOpЂ*sequential_2/dense_8/MatMul/ReadVariableOpЂ,sequential_2/dense_8/MatMul_1/ReadVariableOp
.tf.__operators__.getitem_5/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       
0tf.__operators__.getitem_5/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        
0tf.__operators__.getitem_5/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      а
(tf.__operators__.getitem_5/strided_sliceStridedSliceinputs7tf.__operators__.getitem_5/strided_slice/stack:output:09tf.__operators__.getitem_5/strided_slice/stack_1:output:09tf.__operators__.getitem_5/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*

begin_mask*
end_mask
.tf.__operators__.getitem_4/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        
0tf.__operators__.getitem_4/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       
0tf.__operators__.getitem_4/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      а
(tf.__operators__.getitem_4/strided_sliceStridedSliceinputs7tf.__operators__.getitem_4/strided_slice/stack:output:09tf.__operators__.getitem_4/strided_slice/stack_1:output:09tf.__operators__.getitem_4/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*

begin_mask*
end_mask
*sequential_2/dense_6/MatMul/ReadVariableOpReadVariableOp3sequential_2_dense_6_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0О
sequential_2/dense_6/MatMulMatMul1tf.__operators__.getitem_4/strided_slice:output:02sequential_2/dense_6/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
+sequential_2/dense_6/BiasAdd/ReadVariableOpReadVariableOp4sequential_2_dense_6_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0Е
sequential_2/dense_6/BiasAddBiasAdd%sequential_2/dense_6/MatMul:product:03sequential_2/dense_6/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2z
sequential_2/dense_6/ReluRelu%sequential_2/dense_6/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
*sequential_2/dense_7/MatMul/ReadVariableOpReadVariableOp3sequential_2_dense_7_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Д
sequential_2/dense_7/MatMulMatMul'sequential_2/dense_6/Relu:activations:02sequential_2/dense_7/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
+sequential_2/dense_7/BiasAdd/ReadVariableOpReadVariableOp4sequential_2_dense_7_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0Е
sequential_2/dense_7/BiasAddBiasAdd%sequential_2/dense_7/MatMul:product:03sequential_2/dense_7/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2z
sequential_2/dense_7/ReluRelu%sequential_2/dense_7/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
*sequential_2/dense_8/MatMul/ReadVariableOpReadVariableOp3sequential_2_dense_8_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0Д
sequential_2/dense_8/MatMulMatMul'sequential_2/dense_7/Relu:activations:02sequential_2/dense_8/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
+sequential_2/dense_8/BiasAdd/ReadVariableOpReadVariableOp4sequential_2_dense_8_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0Е
sequential_2/dense_8/BiasAddBiasAdd%sequential_2/dense_8/MatMul:product:03sequential_2/dense_8/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ 
,sequential_2/dense_6/MatMul_1/ReadVariableOpReadVariableOp3sequential_2_dense_6_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0Т
sequential_2/dense_6/MatMul_1MatMul1tf.__operators__.getitem_5/strided_slice:output:04sequential_2/dense_6/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
-sequential_2/dense_6/BiasAdd_1/ReadVariableOpReadVariableOp4sequential_2_dense_6_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0Л
sequential_2/dense_6/BiasAdd_1BiasAdd'sequential_2/dense_6/MatMul_1:product:05sequential_2/dense_6/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2~
sequential_2/dense_6/Relu_1Relu'sequential_2/dense_6/BiasAdd_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
,sequential_2/dense_7/MatMul_1/ReadVariableOpReadVariableOp3sequential_2_dense_7_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0К
sequential_2/dense_7/MatMul_1MatMul)sequential_2/dense_6/Relu_1:activations:04sequential_2/dense_7/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
-sequential_2/dense_7/BiasAdd_1/ReadVariableOpReadVariableOp4sequential_2_dense_7_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0Л
sequential_2/dense_7/BiasAdd_1BiasAdd'sequential_2/dense_7/MatMul_1:product:05sequential_2/dense_7/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2~
sequential_2/dense_7/Relu_1Relu'sequential_2/dense_7/BiasAdd_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
,sequential_2/dense_8/MatMul_1/ReadVariableOpReadVariableOp3sequential_2_dense_8_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0К
sequential_2/dense_8/MatMul_1MatMul)sequential_2/dense_7/Relu_1:activations:04sequential_2/dense_8/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ
-sequential_2/dense_8/BiasAdd_1/ReadVariableOpReadVariableOp4sequential_2_dense_8_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0Л
sequential_2/dense_8/BiasAdd_1BiasAdd'sequential_2/dense_8/MatMul_1:product:05sequential_2/dense_8/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџГ
tf.stack_2/stackPack%sequential_2/dense_8/BiasAdd:output:0'sequential_2/dense_8/BiasAdd_1:output:0*
N*
T0*+
_output_shapes
:џџџџџџџџџ*

axisl
IdentityIdentitytf.stack_2/stack:output:0^NoOp*
T0*+
_output_shapes
:џџџџџџџџџє
NoOpNoOp,^sequential_2/dense_6/BiasAdd/ReadVariableOp.^sequential_2/dense_6/BiasAdd_1/ReadVariableOp+^sequential_2/dense_6/MatMul/ReadVariableOp-^sequential_2/dense_6/MatMul_1/ReadVariableOp,^sequential_2/dense_7/BiasAdd/ReadVariableOp.^sequential_2/dense_7/BiasAdd_1/ReadVariableOp+^sequential_2/dense_7/MatMul/ReadVariableOp-^sequential_2/dense_7/MatMul_1/ReadVariableOp,^sequential_2/dense_8/BiasAdd/ReadVariableOp.^sequential_2/dense_8/BiasAdd_1/ReadVariableOp+^sequential_2/dense_8/MatMul/ReadVariableOp-^sequential_2/dense_8/MatMul_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 2Z
+sequential_2/dense_6/BiasAdd/ReadVariableOp+sequential_2/dense_6/BiasAdd/ReadVariableOp2^
-sequential_2/dense_6/BiasAdd_1/ReadVariableOp-sequential_2/dense_6/BiasAdd_1/ReadVariableOp2X
*sequential_2/dense_6/MatMul/ReadVariableOp*sequential_2/dense_6/MatMul/ReadVariableOp2\
,sequential_2/dense_6/MatMul_1/ReadVariableOp,sequential_2/dense_6/MatMul_1/ReadVariableOp2Z
+sequential_2/dense_7/BiasAdd/ReadVariableOp+sequential_2/dense_7/BiasAdd/ReadVariableOp2^
-sequential_2/dense_7/BiasAdd_1/ReadVariableOp-sequential_2/dense_7/BiasAdd_1/ReadVariableOp2X
*sequential_2/dense_7/MatMul/ReadVariableOp*sequential_2/dense_7/MatMul/ReadVariableOp2\
,sequential_2/dense_7/MatMul_1/ReadVariableOp,sequential_2/dense_7/MatMul_1/ReadVariableOp2Z
+sequential_2/dense_8/BiasAdd/ReadVariableOp+sequential_2/dense_8/BiasAdd/ReadVariableOp2^
-sequential_2/dense_8/BiasAdd_1/ReadVariableOp-sequential_2/dense_8/BiasAdd_1/ReadVariableOp2X
*sequential_2/dense_8/MatMul/ReadVariableOp*sequential_2/dense_8/MatMul/ReadVariableOp2\
,sequential_2/dense_8/MatMul_1/ReadVariableOp,sequential_2/dense_8/MatMul_1/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
Х4
ё	
 __inference__traced_save_1184650
file_prefix-
)savev2_dense_6_kernel_read_readvariableop+
'savev2_dense_6_bias_read_readvariableop-
)savev2_dense_7_kernel_read_readvariableop+
'savev2_dense_7_bias_read_readvariableop-
)savev2_dense_8_kernel_read_readvariableop+
'savev2_dense_8_bias_read_readvariableop(
$savev2_iteration_read_readvariableop	,
(savev2_learning_rate_read_readvariableop4
0savev2_adam_m_dense_6_kernel_read_readvariableop4
0savev2_adam_v_dense_6_kernel_read_readvariableop2
.savev2_adam_m_dense_6_bias_read_readvariableop2
.savev2_adam_v_dense_6_bias_read_readvariableop4
0savev2_adam_m_dense_7_kernel_read_readvariableop4
0savev2_adam_v_dense_7_kernel_read_readvariableop2
.savev2_adam_m_dense_7_bias_read_readvariableop2
.savev2_adam_v_dense_7_bias_read_readvariableop4
0savev2_adam_m_dense_8_kernel_read_readvariableop4
0savev2_adam_v_dense_8_kernel_read_readvariableop2
.savev2_adam_m_dense_8_bias_read_readvariableop2
.savev2_adam_v_dense_8_bias_read_readvariableop&
"savev2_total_1_read_readvariableop&
"savev2_count_1_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop
savev2_const

identity_1ЂMergeV2Checkpointsw
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*Z
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.parta
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B
_temp/part
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: f

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: L

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :f
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : 
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
:  

SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*Щ	
valueП	BМ	B&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUEB3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*E
value<B:B B B B B B B B B B B B B B B B B B B B B B B B B 

SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0)savev2_dense_6_kernel_read_readvariableop'savev2_dense_6_bias_read_readvariableop)savev2_dense_7_kernel_read_readvariableop'savev2_dense_7_bias_read_readvariableop)savev2_dense_8_kernel_read_readvariableop'savev2_dense_8_bias_read_readvariableop$savev2_iteration_read_readvariableop(savev2_learning_rate_read_readvariableop0savev2_adam_m_dense_6_kernel_read_readvariableop0savev2_adam_v_dense_6_kernel_read_readvariableop.savev2_adam_m_dense_6_bias_read_readvariableop.savev2_adam_v_dense_6_bias_read_readvariableop0savev2_adam_m_dense_7_kernel_read_readvariableop0savev2_adam_v_dense_7_kernel_read_readvariableop.savev2_adam_m_dense_7_bias_read_readvariableop.savev2_adam_v_dense_7_bias_read_readvariableop0savev2_adam_m_dense_8_kernel_read_readvariableop0savev2_adam_v_dense_8_kernel_read_readvariableop.savev2_adam_m_dense_8_bias_read_readvariableop.savev2_adam_v_dense_8_bias_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableopsavev2_const"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *'
dtypes
2	
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:Г
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 f
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: Q

Identity_1IdentityIdentity:output:0^NoOp*
T0*
_output_shapes
: [
NoOpNoOp^MergeV2Checkpoints*"
_acd_function_control_output(*
_output_shapes
 "!

identity_1Identity_1:output:0*Е
_input_shapesЃ
 : :2:2:22:2:2:: : :2:2:2:2:22:22:2:2:2:2::: : : : : 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:$ 

_output_shapes

:2: 

_output_shapes
:2:$ 

_output_shapes

:22: 

_output_shapes
:2:$ 

_output_shapes

:2: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :$	 

_output_shapes

:2:$
 

_output_shapes

:2: 

_output_shapes
:2: 

_output_shapes
:2:$ 

_output_shapes

:22:$ 

_output_shapes

:22: 

_output_shapes
:2: 

_output_shapes
:2:$ 

_output_shapes

:2:$ 

_output_shapes

:2: 

_output_shapes
:: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
	

.__inference_sequential_2_layer_call_fn_1183995
dense_6_input
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityЂStatefulPartitionedCall
StatefulPartitionedCallStatefulPartitionedCalldense_6_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *R
fMRK
I__inference_sequential_2_layer_call_and_return_conditional_losses_1183963o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџ`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:V R
'
_output_shapes
:џџџџџџџџџ
'
_user_specified_namedense_6_input


I__inference_sequential_2_layer_call_and_return_conditional_losses_1184033
dense_6_input!
dense_6_1184017:2
dense_6_1184019:2!
dense_7_1184022:22
dense_7_1184024:2!
dense_8_1184027:2
dense_8_1184029:
identityЂdense_6/StatefulPartitionedCallЂdense_7/StatefulPartitionedCallЂdense_8/StatefulPartitionedCallћ
dense_6/StatefulPartitionedCallStatefulPartitionedCalldense_6_inputdense_6_1184017dense_6_1184019*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_dense_6_layer_call_and_return_conditional_losses_1183840
dense_7/StatefulPartitionedCallStatefulPartitionedCall(dense_6/StatefulPartitionedCall:output:0dense_7_1184022dense_7_1184024*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_dense_7_layer_call_and_return_conditional_losses_1183857
dense_8/StatefulPartitionedCallStatefulPartitionedCall(dense_7/StatefulPartitionedCall:output:0dense_8_1184027dense_8_1184029*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_dense_8_layer_call_and_return_conditional_losses_1183873w
IdentityIdentity(dense_8/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџЌ
NoOpNoOp ^dense_6/StatefulPartitionedCall ^dense_7/StatefulPartitionedCall ^dense_8/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 2B
dense_6/StatefulPartitionedCalldense_6/StatefulPartitionedCall2B
dense_7/StatefulPartitionedCalldense_7/StatefulPartitionedCall2B
dense_8/StatefulPartitionedCalldense_8/StatefulPartitionedCall:V R
'
_output_shapes
:џџџџџџџџџ
'
_user_specified_namedense_6_input


ѕ
D__inference_dense_7_layer_call_and_return_conditional_losses_1183857

inputs0
matmul_readvariableop_resource:22-
biasadd_readvariableop_resource:2
identityЂBiasAdd/ReadVariableOpЂMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:22*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџ2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:џџџџџџџџџ2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:џџџџџџџџџ2
 
_user_specified_nameinputs


I__inference_sequential_2_layer_call_and_return_conditional_losses_1183880

inputs!
dense_6_1183841:2
dense_6_1183843:2!
dense_7_1183858:22
dense_7_1183860:2!
dense_8_1183874:2
dense_8_1183876:
identityЂdense_6/StatefulPartitionedCallЂdense_7/StatefulPartitionedCallЂdense_8/StatefulPartitionedCallє
dense_6/StatefulPartitionedCallStatefulPartitionedCallinputsdense_6_1183841dense_6_1183843*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_dense_6_layer_call_and_return_conditional_losses_1183840
dense_7/StatefulPartitionedCallStatefulPartitionedCall(dense_6/StatefulPartitionedCall:output:0dense_7_1183858dense_7_1183860*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_dense_7_layer_call_and_return_conditional_losses_1183857
dense_8/StatefulPartitionedCallStatefulPartitionedCall(dense_7/StatefulPartitionedCall:output:0dense_8_1183874dense_8_1183876*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8 *M
fHRF
D__inference_dense_8_layer_call_and_return_conditional_losses_1183873w
IdentityIdentity(dense_8/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:џџџџџџџџџЌ
NoOpNoOp ^dense_6/StatefulPartitionedCall ^dense_7/StatefulPartitionedCall ^dense_8/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 2B
dense_6/StatefulPartitionedCalldense_6/StatefulPartitionedCall2B
dense_7/StatefulPartitionedCalldense_7/StatefulPartitionedCall2B
dense_8/StatefulPartitionedCalldense_8/StatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
Л

D__inference_model_2_layer_call_and_return_conditional_losses_1184137

inputs&
sequential_2_1184115:2"
sequential_2_1184117:2&
sequential_2_1184119:22"
sequential_2_1184121:2&
sequential_2_1184123:2"
sequential_2_1184125:
identityЂ$sequential_2/StatefulPartitionedCallЂ&sequential_2/StatefulPartitionedCall_1
.tf.__operators__.getitem_5/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       
0tf.__operators__.getitem_5/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        
0tf.__operators__.getitem_5/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      а
(tf.__operators__.getitem_5/strided_sliceStridedSliceinputs7tf.__operators__.getitem_5/strided_slice/stack:output:09tf.__operators__.getitem_5/strided_slice/stack_1:output:09tf.__operators__.getitem_5/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*

begin_mask*
end_mask
.tf.__operators__.getitem_4/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        
0tf.__operators__.getitem_4/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       
0tf.__operators__.getitem_4/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      а
(tf.__operators__.getitem_4/strided_sliceStridedSliceinputs7tf.__operators__.getitem_4/strided_slice/stack:output:09tf.__operators__.getitem_4/strided_slice/stack_1:output:09tf.__operators__.getitem_4/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*

begin_mask*
end_mask
$sequential_2/StatefulPartitionedCallStatefulPartitionedCall1tf.__operators__.getitem_4/strided_slice:output:0sequential_2_1184115sequential_2_1184117sequential_2_1184119sequential_2_1184121sequential_2_1184123sequential_2_1184125*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *R
fMRK
I__inference_sequential_2_layer_call_and_return_conditional_losses_1183963
&sequential_2/StatefulPartitionedCall_1StatefulPartitionedCall1tf.__operators__.getitem_5/strided_slice:output:0sequential_2_1184115sequential_2_1184117sequential_2_1184119sequential_2_1184121sequential_2_1184123sequential_2_1184125*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8 *R
fMRK
I__inference_sequential_2_layer_call_and_return_conditional_losses_1183963У
tf.stack_2/stackPack-sequential_2/StatefulPartitionedCall:output:0/sequential_2/StatefulPartitionedCall_1:output:0*
N*
T0*+
_output_shapes
:џџџџџџџџџ*

axisl
IdentityIdentitytf.stack_2/stack:output:0^NoOp*
T0*+
_output_shapes
:џџџџџџџџџ
NoOpNoOp%^sequential_2/StatefulPartitionedCall'^sequential_2/StatefulPartitionedCall_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:џџџџџџџџџ: : : : : : 2L
$sequential_2/StatefulPartitionedCall$sequential_2/StatefulPartitionedCall2P
&sequential_2/StatefulPartitionedCall_1&sequential_2/StatefulPartitionedCall_1:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
­
L
$__inference__update_step_xla_1184394
gradient
variable:2*
_XlaMustCompile(*(
_construction_contextkEagerRuntime*
_input_shapes

:2: *
	_noinline(:D @

_output_shapes
:2
"
_user_specified_name
gradient:($
"
_user_specified_name
variable
­
L
$__inference__update_step_xla_1184404
gradient
variable:2*
_XlaMustCompile(*(
_construction_contextkEagerRuntime*
_input_shapes

:2: *
	_noinline(:D @

_output_shapes
:2
"
_user_specified_name
gradient:($
"
_user_specified_name
variable"
L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*Б
serving_default
;
input_30
serving_default_input_3:0џџџџџџџџџB

tf.stack_24
StatefulPartitionedCall:0џџџџџџџџџtensorflow/serving/predict:АЕ
О
layer-0
layer-1
layer-2
layer_with_weights-0
layer-3
layer-4
	variables
trainable_variables
regularization_losses
		keras_api

__call__
*&call_and_return_all_conditional_losses
_default_save_signature
	optimizer

signatures"
_tf_keras_network
"
_tf_keras_input_layer
(
	keras_api"
_tf_keras_layer
(
	keras_api"
_tf_keras_layer

layer_with_weights-0
layer-0
layer_with_weights-1
layer-1
layer_with_weights-2
layer-2
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses"
_tf_keras_sequential
(
	keras_api"
_tf_keras_layer
J
0
1
2
3
4
 5"
trackable_list_wrapper
J
0
1
2
3
4
 5"
trackable_list_wrapper
 "
trackable_list_wrapper
Ъ
!non_trainable_variables

"layers
#metrics
$layer_regularization_losses
%layer_metrics
	variables
trainable_variables
regularization_losses

__call__
_default_save_signature
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses"
_generic_user_object
й
&trace_0
'trace_1
(trace_2
)trace_32ю
)__inference_model_2_layer_call_fn_1184085
)__inference_model_2_layer_call_fn_1184273
)__inference_model_2_layer_call_fn_1184290
)__inference_model_2_layer_call_fn_1184169П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 z&trace_0z'trace_1z(trace_2z)trace_3
Х
*trace_0
+trace_1
,trace_2
-trace_32к
D__inference_model_2_layer_call_and_return_conditional_losses_1184337
D__inference_model_2_layer_call_and_return_conditional_losses_1184384
D__inference_model_2_layer_call_and_return_conditional_losses_1184202
D__inference_model_2_layer_call_and_return_conditional_losses_1184235П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 z*trace_0z+trace_1z,trace_2z-trace_3
ЭBЪ
"__inference__wrapped_model_1183822input_3"
В
FullArgSpec
args 
varargsjargs
varkwjkwargs
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 

.
_variables
/_iterations
0_learning_rate
1_index_dict
2
_momentums
3_velocities
4_update_step_xla"
experimentalOptimizer
,
5serving_default"
signature_map
"
_generic_user_object
"
_generic_user_object
Л
6	variables
7trainable_variables
8regularization_losses
9	keras_api
:__call__
*;&call_and_return_all_conditional_losses

kernel
bias"
_tf_keras_layer
Л
<	variables
=trainable_variables
>regularization_losses
?	keras_api
@__call__
*A&call_and_return_all_conditional_losses

kernel
bias"
_tf_keras_layer
Л
B	variables
Ctrainable_variables
Dregularization_losses
E	keras_api
F__call__
*G&call_and_return_all_conditional_losses

kernel
 bias"
_tf_keras_layer
J
0
1
2
3
4
 5"
trackable_list_wrapper
J
0
1
2
3
4
 5"
trackable_list_wrapper
 "
trackable_list_wrapper
­
Hnon_trainable_variables

Ilayers
Jmetrics
Klayer_regularization_losses
Llayer_metrics
	variables
trainable_variables
regularization_losses
__call__
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses"
_generic_user_object
э
Mtrace_0
Ntrace_1
Otrace_2
Ptrace_32
.__inference_sequential_2_layer_call_fn_1183895
.__inference_sequential_2_layer_call_fn_1184431
.__inference_sequential_2_layer_call_fn_1184448
.__inference_sequential_2_layer_call_fn_1183995П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 zMtrace_0zNtrace_1zOtrace_2zPtrace_3
й
Qtrace_0
Rtrace_1
Strace_2
Ttrace_32ю
I__inference_sequential_2_layer_call_and_return_conditional_losses_1184472
I__inference_sequential_2_layer_call_and_return_conditional_losses_1184496
I__inference_sequential_2_layer_call_and_return_conditional_losses_1184014
I__inference_sequential_2_layer_call_and_return_conditional_losses_1184033П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 zQtrace_0zRtrace_1zStrace_2zTtrace_3
"
_generic_user_object
 :22dense_6/kernel
:22dense_6/bias
 :222dense_7/kernel
:22dense_7/bias
 :22dense_8/kernel
:2dense_8/bias
 "
trackable_list_wrapper
C
0
1
2
3
4"
trackable_list_wrapper
.
U0
V1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
ћBј
)__inference_model_2_layer_call_fn_1184085input_3"П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
њBї
)__inference_model_2_layer_call_fn_1184273inputs"П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
њBї
)__inference_model_2_layer_call_fn_1184290inputs"П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
ћBј
)__inference_model_2_layer_call_fn_1184169input_3"П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
B
D__inference_model_2_layer_call_and_return_conditional_losses_1184337inputs"П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
B
D__inference_model_2_layer_call_and_return_conditional_losses_1184384inputs"П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
B
D__inference_model_2_layer_call_and_return_conditional_losses_1184202input_3"П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
B
D__inference_model_2_layer_call_and_return_conditional_losses_1184235input_3"П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
~
/0
W1
X2
Y3
Z4
[5
\6
]7
^8
_9
`10
a11
b12"
trackable_list_wrapper
:	 2	iteration
: 2learning_rate
 "
trackable_dict_wrapper
J
W0
Y1
[2
]3
_4
a5"
trackable_list_wrapper
J
X0
Z1
\2
^3
`4
b5"
trackable_list_wrapper
П
ctrace_0
dtrace_1
etrace_2
ftrace_3
gtrace_4
htrace_52 
$__inference__update_step_xla_1184389
$__inference__update_step_xla_1184394
$__inference__update_step_xla_1184399
$__inference__update_step_xla_1184404
$__inference__update_step_xla_1184409
$__inference__update_step_xla_1184414Й
ЎВЊ
FullArgSpec2
args*'
jself

jgradient

jvariable
jkey
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 0zctrace_0zdtrace_1zetrace_2zftrace_3zgtrace_4zhtrace_5
ЬBЩ
%__inference_signature_wrapper_1184256input_3"
В
FullArgSpec
args 
varargs
 
varkwjkwargs
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
­
inon_trainable_variables

jlayers
kmetrics
llayer_regularization_losses
mlayer_metrics
6	variables
7trainable_variables
8regularization_losses
:__call__
*;&call_and_return_all_conditional_losses
&;"call_and_return_conditional_losses"
_generic_user_object
э
ntrace_02а
)__inference_dense_6_layer_call_fn_1184505Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 zntrace_0

otrace_02ы
D__inference_dense_6_layer_call_and_return_conditional_losses_1184516Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 zotrace_0
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
­
pnon_trainable_variables

qlayers
rmetrics
slayer_regularization_losses
tlayer_metrics
<	variables
=trainable_variables
>regularization_losses
@__call__
*A&call_and_return_all_conditional_losses
&A"call_and_return_conditional_losses"
_generic_user_object
э
utrace_02а
)__inference_dense_7_layer_call_fn_1184525Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 zutrace_0

vtrace_02ы
D__inference_dense_7_layer_call_and_return_conditional_losses_1184536Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 zvtrace_0
.
0
 1"
trackable_list_wrapper
.
0
 1"
trackable_list_wrapper
 "
trackable_list_wrapper
­
wnon_trainable_variables

xlayers
ymetrics
zlayer_regularization_losses
{layer_metrics
B	variables
Ctrainable_variables
Dregularization_losses
F__call__
*G&call_and_return_all_conditional_losses
&G"call_and_return_conditional_losses"
_generic_user_object
э
|trace_02а
)__inference_dense_8_layer_call_fn_1184545Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 z|trace_0

}trace_02ы
D__inference_dense_8_layer_call_and_return_conditional_losses_1184555Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 z}trace_0
 "
trackable_list_wrapper
5
0
1
2"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
B
.__inference_sequential_2_layer_call_fn_1183895dense_6_input"П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
џBќ
.__inference_sequential_2_layer_call_fn_1184431inputs"П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
џBќ
.__inference_sequential_2_layer_call_fn_1184448inputs"П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
B
.__inference_sequential_2_layer_call_fn_1183995dense_6_input"П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
B
I__inference_sequential_2_layer_call_and_return_conditional_losses_1184472inputs"П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
B
I__inference_sequential_2_layer_call_and_return_conditional_losses_1184496inputs"П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
ЁB
I__inference_sequential_2_layer_call_and_return_conditional_losses_1184014dense_6_input"П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
ЁB
I__inference_sequential_2_layer_call_and_return_conditional_losses_1184033dense_6_input"П
ЖВВ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
P
~	variables
	keras_api

total

count"
_tf_keras_metric
c
	variables
	keras_api

total

count

_fn_kwargs"
_tf_keras_metric
%:#22Adam/m/dense_6/kernel
%:#22Adam/v/dense_6/kernel
:22Adam/m/dense_6/bias
:22Adam/v/dense_6/bias
%:#222Adam/m/dense_7/kernel
%:#222Adam/v/dense_7/kernel
:22Adam/m/dense_7/bias
:22Adam/v/dense_7/bias
%:#22Adam/m/dense_8/kernel
%:#22Adam/v/dense_8/kernel
:2Adam/m/dense_8/bias
:2Adam/v/dense_8/bias
љBі
$__inference__update_step_xla_1184389gradientvariable"З
ЎВЊ
FullArgSpec2
args*'
jself

jgradient

jvariable
jkey
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
љBі
$__inference__update_step_xla_1184394gradientvariable"З
ЎВЊ
FullArgSpec2
args*'
jself

jgradient

jvariable
jkey
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
љBі
$__inference__update_step_xla_1184399gradientvariable"З
ЎВЊ
FullArgSpec2
args*'
jself

jgradient

jvariable
jkey
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
љBі
$__inference__update_step_xla_1184404gradientvariable"З
ЎВЊ
FullArgSpec2
args*'
jself

jgradient

jvariable
jkey
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
љBі
$__inference__update_step_xla_1184409gradientvariable"З
ЎВЊ
FullArgSpec2
args*'
jself

jgradient

jvariable
jkey
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
љBі
$__inference__update_step_xla_1184414gradientvariable"З
ЎВЊ
FullArgSpec2
args*'
jself

jgradient

jvariable
jkey
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
нBк
)__inference_dense_6_layer_call_fn_1184505inputs"Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
јBѕ
D__inference_dense_6_layer_call_and_return_conditional_losses_1184516inputs"Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
нBк
)__inference_dense_7_layer_call_fn_1184525inputs"Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
јBѕ
D__inference_dense_7_layer_call_and_return_conditional_losses_1184536inputs"Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
нBк
)__inference_dense_8_layer_call_fn_1184545inputs"Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
јBѕ
D__inference_dense_8_layer_call_and_return_conditional_losses_1184555inputs"Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
0
0
1"
trackable_list_wrapper
-
~	variables"
_generic_user_object
:  (2total
:  (2count
0
0
1"
trackable_list_wrapper
.
	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
$__inference__update_step_xla_1184389nhЂe
^Ђ[

gradient2
41	Ђ
њ2

p
` VariableSpec 
`рЗкџцп?
Њ "
 
$__inference__update_step_xla_1184394f`Ђ]
VЂS

gradient2
0-	Ђ
њ2

p
` VariableSpec 
` лџцп?
Њ "
 
$__inference__update_step_xla_1184399nhЂe
^Ђ[

gradient22
41	Ђ
њ22

p
` VariableSpec 
` Ікџцп?
Њ "
 
$__inference__update_step_xla_1184404f`Ђ]
VЂS

gradient2
0-	Ђ
њ2

p
` VariableSpec 
`Њкџцп?
Њ "
 
$__inference__update_step_xla_1184409nhЂe
^Ђ[

gradient2
41	Ђ
њ2

p
` VariableSpec 
`рѓкџцп?
Њ "
 
$__inference__update_step_xla_1184414f`Ђ]
VЂS

gradient
0-	Ђ
њ

p
` VariableSpec 
`длџцп?
Њ "
 
"__inference__wrapped_model_1183822w 0Ђ-
&Ђ#
!
input_3џџџџџџџџџ
Њ ";Њ8
6

tf.stack_2(%

tf_stack_2џџџџџџџџџЋ
D__inference_dense_6_layer_call_and_return_conditional_losses_1184516c/Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ ",Ђ)
"
tensor_0џџџџџџџџџ2
 
)__inference_dense_6_layer_call_fn_1184505X/Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "!
unknownџџџџџџџџџ2Ћ
D__inference_dense_7_layer_call_and_return_conditional_losses_1184536c/Ђ,
%Ђ"
 
inputsџџџџџџџџџ2
Њ ",Ђ)
"
tensor_0џџџџџџџџџ2
 
)__inference_dense_7_layer_call_fn_1184525X/Ђ,
%Ђ"
 
inputsџџџџџџџџџ2
Њ "!
unknownџџџџџџџџџ2Ћ
D__inference_dense_8_layer_call_and_return_conditional_losses_1184555c /Ђ,
%Ђ"
 
inputsџџџџџџџџџ2
Њ ",Ђ)
"
tensor_0џџџџџџџџџ
 
)__inference_dense_8_layer_call_fn_1184545X /Ђ,
%Ђ"
 
inputsџџџџџџџџџ2
Њ "!
unknownџџџџџџџџџМ
D__inference_model_2_layer_call_and_return_conditional_losses_1184202t 8Ђ5
.Ђ+
!
input_3џџџџџџџџџ
p 

 
Њ "0Ђ-
&#
tensor_0џџџџџџџџџ
 М
D__inference_model_2_layer_call_and_return_conditional_losses_1184235t 8Ђ5
.Ђ+
!
input_3џџџџџџџџџ
p

 
Њ "0Ђ-
&#
tensor_0џџџџџџџџџ
 Л
D__inference_model_2_layer_call_and_return_conditional_losses_1184337s 7Ђ4
-Ђ*
 
inputsџџџџџџџџџ
p 

 
Њ "0Ђ-
&#
tensor_0џџџџџџџџџ
 Л
D__inference_model_2_layer_call_and_return_conditional_losses_1184384s 7Ђ4
-Ђ*
 
inputsџџџџџџџџџ
p

 
Њ "0Ђ-
&#
tensor_0џџџџџџџџџ
 
)__inference_model_2_layer_call_fn_1184085i 8Ђ5
.Ђ+
!
input_3џџџџџџџџџ
p 

 
Њ "%"
unknownџџџџџџџџџ
)__inference_model_2_layer_call_fn_1184169i 8Ђ5
.Ђ+
!
input_3џџџџџџџџџ
p

 
Њ "%"
unknownџџџџџџџџџ
)__inference_model_2_layer_call_fn_1184273h 7Ђ4
-Ђ*
 
inputsџџџџџџџџџ
p 

 
Њ "%"
unknownџџџџџџџџџ
)__inference_model_2_layer_call_fn_1184290h 7Ђ4
-Ђ*
 
inputsџџџџџџџџџ
p

 
Њ "%"
unknownџџџџџџџџџУ
I__inference_sequential_2_layer_call_and_return_conditional_losses_1184014v >Ђ;
4Ђ1
'$
dense_6_inputџџџџџџџџџ
p 

 
Њ ",Ђ)
"
tensor_0џџџџџџџџџ
 У
I__inference_sequential_2_layer_call_and_return_conditional_losses_1184033v >Ђ;
4Ђ1
'$
dense_6_inputџџџџџџџџџ
p

 
Њ ",Ђ)
"
tensor_0џџџџџџџџџ
 М
I__inference_sequential_2_layer_call_and_return_conditional_losses_1184472o 7Ђ4
-Ђ*
 
inputsџџџџџџџџџ
p 

 
Њ ",Ђ)
"
tensor_0џџџџџџџџџ
 М
I__inference_sequential_2_layer_call_and_return_conditional_losses_1184496o 7Ђ4
-Ђ*
 
inputsџџџџџџџџџ
p

 
Њ ",Ђ)
"
tensor_0џџџџџџџџџ
 
.__inference_sequential_2_layer_call_fn_1183895k >Ђ;
4Ђ1
'$
dense_6_inputџџџџџџџџџ
p 

 
Њ "!
unknownџџџџџџџџџ
.__inference_sequential_2_layer_call_fn_1183995k >Ђ;
4Ђ1
'$
dense_6_inputџџџџџџџџџ
p

 
Њ "!
unknownџџџџџџџџџ
.__inference_sequential_2_layer_call_fn_1184431d 7Ђ4
-Ђ*
 
inputsџџџџџџџџџ
p 

 
Њ "!
unknownџџџџџџџџџ
.__inference_sequential_2_layer_call_fn_1184448d 7Ђ4
-Ђ*
 
inputsџџџџџџџџџ
p

 
Њ "!
unknownџџџџџџџџџЌ
%__inference_signature_wrapper_1184256 ;Ђ8
Ђ 
1Њ.
,
input_3!
input_3џџџџџџџџџ";Њ8
6

tf.stack_2(%

tf_stack_2џџџџџџџџџ