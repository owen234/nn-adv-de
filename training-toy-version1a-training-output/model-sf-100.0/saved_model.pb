▄С
█Й
^
AssignVariableOp
resource
value"dtype"
dtypetype"
validate_shapebool( ѕ
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
є
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool("
allow_missing_filesbool( ѕ
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
dtypetypeѕ
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
list(type)(0ѕ
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0ѕ
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
┴
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
executor_typestring ѕе
@
StaticRegexFullMatch	
input

output
"
patternstring
э
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
ќ
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 ѕ"serve*2.11.02unknown8БЮ
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
ђ
Adam/v/dense_11/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/v/dense_11/bias
y
(Adam/v/dense_11/bias/Read/ReadVariableOpReadVariableOpAdam/v/dense_11/bias*
_output_shapes
:*
dtype0
ђ
Adam/m/dense_11/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/m/dense_11/bias
y
(Adam/m/dense_11/bias/Read/ReadVariableOpReadVariableOpAdam/m/dense_11/bias*
_output_shapes
:*
dtype0
ѕ
Adam/v/dense_11/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*'
shared_nameAdam/v/dense_11/kernel
Ђ
*Adam/v/dense_11/kernel/Read/ReadVariableOpReadVariableOpAdam/v/dense_11/kernel*
_output_shapes

:2*
dtype0
ѕ
Adam/m/dense_11/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*'
shared_nameAdam/m/dense_11/kernel
Ђ
*Adam/m/dense_11/kernel/Read/ReadVariableOpReadVariableOpAdam/m/dense_11/kernel*
_output_shapes

:2*
dtype0
ђ
Adam/v/dense_10/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*%
shared_nameAdam/v/dense_10/bias
y
(Adam/v/dense_10/bias/Read/ReadVariableOpReadVariableOpAdam/v/dense_10/bias*
_output_shapes
:2*
dtype0
ђ
Adam/m/dense_10/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*%
shared_nameAdam/m/dense_10/bias
y
(Adam/m/dense_10/bias/Read/ReadVariableOpReadVariableOpAdam/m/dense_10/bias*
_output_shapes
:2*
dtype0
ѕ
Adam/v/dense_10/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*'
shared_nameAdam/v/dense_10/kernel
Ђ
*Adam/v/dense_10/kernel/Read/ReadVariableOpReadVariableOpAdam/v/dense_10/kernel*
_output_shapes

:22*
dtype0
ѕ
Adam/m/dense_10/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*'
shared_nameAdam/m/dense_10/kernel
Ђ
*Adam/m/dense_10/kernel/Read/ReadVariableOpReadVariableOpAdam/m/dense_10/kernel*
_output_shapes

:22*
dtype0
~
Adam/v/dense_9/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*$
shared_nameAdam/v/dense_9/bias
w
'Adam/v/dense_9/bias/Read/ReadVariableOpReadVariableOpAdam/v/dense_9/bias*
_output_shapes
:2*
dtype0
~
Adam/m/dense_9/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*$
shared_nameAdam/m/dense_9/bias
w
'Adam/m/dense_9/bias/Read/ReadVariableOpReadVariableOpAdam/m/dense_9/bias*
_output_shapes
:2*
dtype0
є
Adam/v/dense_9/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*&
shared_nameAdam/v/dense_9/kernel

)Adam/v/dense_9/kernel/Read/ReadVariableOpReadVariableOpAdam/v/dense_9/kernel*
_output_shapes

:2*
dtype0
є
Adam/m/dense_9/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*&
shared_nameAdam/m/dense_9/kernel

)Adam/m/dense_9/kernel/Read/ReadVariableOpReadVariableOpAdam/m/dense_9/kernel*
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
r
dense_11/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_namedense_11/bias
k
!dense_11/bias/Read/ReadVariableOpReadVariableOpdense_11/bias*
_output_shapes
:*
dtype0
z
dense_11/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2* 
shared_namedense_11/kernel
s
#dense_11/kernel/Read/ReadVariableOpReadVariableOpdense_11/kernel*
_output_shapes

:2*
dtype0
r
dense_10/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*
shared_namedense_10/bias
k
!dense_10/bias/Read/ReadVariableOpReadVariableOpdense_10/bias*
_output_shapes
:2*
dtype0
z
dense_10/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22* 
shared_namedense_10/kernel
s
#dense_10/kernel/Read/ReadVariableOpReadVariableOpdense_10/kernel*
_output_shapes

:22*
dtype0
p
dense_9/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*
shared_namedense_9/bias
i
 dense_9/bias/Read/ReadVariableOpReadVariableOpdense_9/bias*
_output_shapes
:2*
dtype0
x
dense_9/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*
shared_namedense_9/kernel
q
"dense_9/kernel/Read/ReadVariableOpReadVariableOpdense_9/kernel*
_output_shapes

:2*
dtype0
z
serving_default_input_4Placeholder*'
_output_shapes
:         *
dtype0*
shape:         
Ц
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_4dense_9/kerneldense_9/biasdense_10/kerneldense_10/biasdense_11/kerneldense_11/bias*
Tin
	2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *-
f(R&
$__inference_signature_wrapper_172792

NoOpNoOp
»2
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*Ж1
valueЯ1BП1 Bо1
Д
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
Ё
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
░
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
Ђ
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
д
6	variables
7trainable_variables
8regularization_losses
9	keras_api
:__call__
*;&call_and_return_all_conditional_losses

kernel
bias*
д
<	variables
=trainable_variables
>regularization_losses
?	keras_api
@__call__
*A&call_and_return_all_conditional_losses

kernel
bias*
д
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
Њ
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
VARIABLE_VALUEdense_9/kernel&variables/0/.ATTRIBUTES/VARIABLE_VALUE*
LF
VARIABLE_VALUEdense_9/bias&variables/1/.ATTRIBUTES/VARIABLE_VALUE*
OI
VARIABLE_VALUEdense_10/kernel&variables/2/.ATTRIBUTES/VARIABLE_VALUE*
MG
VARIABLE_VALUEdense_10/bias&variables/3/.ATTRIBUTES/VARIABLE_VALUE*
OI
VARIABLE_VALUEdense_11/kernel&variables/4/.ATTRIBUTES/VARIABLE_VALUE*
MG
VARIABLE_VALUEdense_11/bias&variables/5/.ATTRIBUTES/VARIABLE_VALUE*
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
Њ
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
Њ
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
Њ
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

ђtotal

Ђcount*
M
ѓ	variables
Ѓ	keras_api

ёtotal

Ёcount
є
_fn_kwargs*
`Z
VARIABLE_VALUEAdam/m/dense_9/kernel1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUEAdam/v/dense_9/kernel1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEAdam/m/dense_9/bias1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEAdam/v/dense_9/bias1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/m/dense_10/kernel1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/v/dense_10/kernel1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/m/dense_10/bias1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/v/dense_10/bias1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/m/dense_11/kernel1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUE*
b\
VARIABLE_VALUEAdam/v/dense_11/kernel2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUEAdam/m/dense_11/bias2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUEAdam/v/dense_11/bias2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUE*
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
ђ0
Ђ1*

~	variables*
UO
VARIABLE_VALUEtotal_14keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEcount_14keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE*

ё0
Ё1*

ѓ	variables*
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
└	
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename"dense_9/kernel/Read/ReadVariableOp dense_9/bias/Read/ReadVariableOp#dense_10/kernel/Read/ReadVariableOp!dense_10/bias/Read/ReadVariableOp#dense_11/kernel/Read/ReadVariableOp!dense_11/bias/Read/ReadVariableOpiteration/Read/ReadVariableOp!learning_rate/Read/ReadVariableOp)Adam/m/dense_9/kernel/Read/ReadVariableOp)Adam/v/dense_9/kernel/Read/ReadVariableOp'Adam/m/dense_9/bias/Read/ReadVariableOp'Adam/v/dense_9/bias/Read/ReadVariableOp*Adam/m/dense_10/kernel/Read/ReadVariableOp*Adam/v/dense_10/kernel/Read/ReadVariableOp(Adam/m/dense_10/bias/Read/ReadVariableOp(Adam/v/dense_10/bias/Read/ReadVariableOp*Adam/m/dense_11/kernel/Read/ReadVariableOp*Adam/v/dense_11/kernel/Read/ReadVariableOp(Adam/m/dense_11/bias/Read/ReadVariableOp(Adam/v/dense_11/bias/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOpConst*%
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
GPU2 *0J 8ѓ *(
f#R!
__inference__traced_save_173186
█
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenamedense_9/kerneldense_9/biasdense_10/kerneldense_10/biasdense_11/kerneldense_11/bias	iterationlearning_rateAdam/m/dense_9/kernelAdam/v/dense_9/kernelAdam/m/dense_9/biasAdam/v/dense_9/biasAdam/m/dense_10/kernelAdam/v/dense_10/kernelAdam/m/dense_10/biasAdam/v/dense_10/biasAdam/m/dense_11/kernelAdam/v/dense_11/kernelAdam/m/dense_11/biasAdam/v/dense_11/biastotal_1count_1totalcount*$
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
GPU2 *0J 8ѓ *+
f&R$
"__inference__traced_restore_173268ёф
Е
Ё
H__inference_sequential_3_layer_call_and_return_conditional_losses_173032

inputs8
&dense_9_matmul_readvariableop_resource:25
'dense_9_biasadd_readvariableop_resource:29
'dense_10_matmul_readvariableop_resource:226
(dense_10_biasadd_readvariableop_resource:29
'dense_11_matmul_readvariableop_resource:26
(dense_11_biasadd_readvariableop_resource:
identityѕбdense_10/BiasAdd/ReadVariableOpбdense_10/MatMul/ReadVariableOpбdense_11/BiasAdd/ReadVariableOpбdense_11/MatMul/ReadVariableOpбdense_9/BiasAdd/ReadVariableOpбdense_9/MatMul/ReadVariableOpё
dense_9/MatMul/ReadVariableOpReadVariableOp&dense_9_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0y
dense_9/MatMulMatMulinputs%dense_9/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2ѓ
dense_9/BiasAdd/ReadVariableOpReadVariableOp'dense_9_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0ј
dense_9/BiasAddBiasAdddense_9/MatMul:product:0&dense_9/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2`
dense_9/ReluReludense_9/BiasAdd:output:0*
T0*'
_output_shapes
:         2є
dense_10/MatMul/ReadVariableOpReadVariableOp'dense_10_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Ј
dense_10/MatMulMatMuldense_9/Relu:activations:0&dense_10/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2ё
dense_10/BiasAdd/ReadVariableOpReadVariableOp(dense_10_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0Љ
dense_10/BiasAddBiasAdddense_10/MatMul:product:0'dense_10/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2b
dense_10/ReluReludense_10/BiasAdd:output:0*
T0*'
_output_shapes
:         2є
dense_11/MatMul/ReadVariableOpReadVariableOp'dense_11_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0љ
dense_11/MatMulMatMuldense_10/Relu:activations:0&dense_11/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         ё
dense_11/BiasAdd/ReadVariableOpReadVariableOp(dense_11_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0Љ
dense_11/BiasAddBiasAdddense_11/MatMul:product:0'dense_11/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         h
IdentityIdentitydense_11/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:         Ї
NoOpNoOp ^dense_10/BiasAdd/ReadVariableOp^dense_10/MatMul/ReadVariableOp ^dense_11/BiasAdd/ReadVariableOp^dense_11/MatMul/ReadVariableOp^dense_9/BiasAdd/ReadVariableOp^dense_9/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2B
dense_10/BiasAdd/ReadVariableOpdense_10/BiasAdd/ReadVariableOp2@
dense_10/MatMul/ReadVariableOpdense_10/MatMul/ReadVariableOp2B
dense_11/BiasAdd/ReadVariableOpdense_11/BiasAdd/ReadVariableOp2@
dense_11/MatMul/ReadVariableOpdense_11/MatMul/ReadVariableOp2@
dense_9/BiasAdd/ReadVariableOpdense_9/BiasAdd/ReadVariableOp2>
dense_9/MatMul/ReadVariableOpdense_9/MatMul/ReadVariableOp:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
г
K
#__inference__update_step_xla_172950
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
И
O
#__inference__update_step_xla_172945
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
џ

З
C__inference_dense_9_layer_call_and_return_conditional_losses_172376

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:2
identityѕбBiasAdd/ReadVariableOpбMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:2*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:         2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:         2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:         : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
К
ќ
)__inference_dense_10_layer_call_fn_173061

inputs
unknown:22
	unknown_0:2
identityѕбStatefulPartitionedCallя
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8ѓ *M
fHRF
D__inference_dense_10_layer_call_and_return_conditional_losses_172393o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         2`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:         2: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         2
 
_user_specified_nameinputs
оJ
й
C__inference_model_3_layer_call_and_return_conditional_losses_172873

inputsE
3sequential_3_dense_9_matmul_readvariableop_resource:2B
4sequential_3_dense_9_biasadd_readvariableop_resource:2F
4sequential_3_dense_10_matmul_readvariableop_resource:22C
5sequential_3_dense_10_biasadd_readvariableop_resource:2F
4sequential_3_dense_11_matmul_readvariableop_resource:2C
5sequential_3_dense_11_biasadd_readvariableop_resource:
identityѕб,sequential_3/dense_10/BiasAdd/ReadVariableOpб.sequential_3/dense_10/BiasAdd_1/ReadVariableOpб+sequential_3/dense_10/MatMul/ReadVariableOpб-sequential_3/dense_10/MatMul_1/ReadVariableOpб,sequential_3/dense_11/BiasAdd/ReadVariableOpб.sequential_3/dense_11/BiasAdd_1/ReadVariableOpб+sequential_3/dense_11/MatMul/ReadVariableOpб-sequential_3/dense_11/MatMul_1/ReadVariableOpб+sequential_3/dense_9/BiasAdd/ReadVariableOpб-sequential_3/dense_9/BiasAdd_1/ReadVariableOpб*sequential_3/dense_9/MatMul/ReadVariableOpб,sequential_3/dense_9/MatMul_1/ReadVariableOp
.tf.__operators__.getitem_7/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       Ђ
0tf.__operators__.getitem_7/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        Ђ
0tf.__operators__.getitem_7/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      л
(tf.__operators__.getitem_7/strided_sliceStridedSliceinputs7tf.__operators__.getitem_7/strided_slice/stack:output:09tf.__operators__.getitem_7/strided_slice/stack_1:output:09tf.__operators__.getitem_7/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask
.tf.__operators__.getitem_6/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        Ђ
0tf.__operators__.getitem_6/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       Ђ
0tf.__operators__.getitem_6/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      л
(tf.__operators__.getitem_6/strided_sliceStridedSliceinputs7tf.__operators__.getitem_6/strided_slice/stack:output:09tf.__operators__.getitem_6/strided_slice/stack_1:output:09tf.__operators__.getitem_6/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_maskъ
*sequential_3/dense_9/MatMul/ReadVariableOpReadVariableOp3sequential_3_dense_9_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0Й
sequential_3/dense_9/MatMulMatMul1tf.__operators__.getitem_6/strided_slice:output:02sequential_3/dense_9/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2ю
+sequential_3/dense_9/BiasAdd/ReadVariableOpReadVariableOp4sequential_3_dense_9_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0х
sequential_3/dense_9/BiasAddBiasAdd%sequential_3/dense_9/MatMul:product:03sequential_3/dense_9/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2z
sequential_3/dense_9/ReluRelu%sequential_3/dense_9/BiasAdd:output:0*
T0*'
_output_shapes
:         2а
+sequential_3/dense_10/MatMul/ReadVariableOpReadVariableOp4sequential_3_dense_10_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Х
sequential_3/dense_10/MatMulMatMul'sequential_3/dense_9/Relu:activations:03sequential_3/dense_10/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2ъ
,sequential_3/dense_10/BiasAdd/ReadVariableOpReadVariableOp5sequential_3_dense_10_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0И
sequential_3/dense_10/BiasAddBiasAdd&sequential_3/dense_10/MatMul:product:04sequential_3/dense_10/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2|
sequential_3/dense_10/ReluRelu&sequential_3/dense_10/BiasAdd:output:0*
T0*'
_output_shapes
:         2а
+sequential_3/dense_11/MatMul/ReadVariableOpReadVariableOp4sequential_3_dense_11_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0и
sequential_3/dense_11/MatMulMatMul(sequential_3/dense_10/Relu:activations:03sequential_3/dense_11/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         ъ
,sequential_3/dense_11/BiasAdd/ReadVariableOpReadVariableOp5sequential_3_dense_11_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0И
sequential_3/dense_11/BiasAddBiasAdd&sequential_3/dense_11/MatMul:product:04sequential_3/dense_11/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         а
,sequential_3/dense_9/MatMul_1/ReadVariableOpReadVariableOp3sequential_3_dense_9_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0┬
sequential_3/dense_9/MatMul_1MatMul1tf.__operators__.getitem_7/strided_slice:output:04sequential_3/dense_9/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2ъ
-sequential_3/dense_9/BiasAdd_1/ReadVariableOpReadVariableOp4sequential_3_dense_9_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0╗
sequential_3/dense_9/BiasAdd_1BiasAdd'sequential_3/dense_9/MatMul_1:product:05sequential_3/dense_9/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2~
sequential_3/dense_9/Relu_1Relu'sequential_3/dense_9/BiasAdd_1:output:0*
T0*'
_output_shapes
:         2б
-sequential_3/dense_10/MatMul_1/ReadVariableOpReadVariableOp4sequential_3_dense_10_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0╝
sequential_3/dense_10/MatMul_1MatMul)sequential_3/dense_9/Relu_1:activations:05sequential_3/dense_10/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2а
.sequential_3/dense_10/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_3_dense_10_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0Й
sequential_3/dense_10/BiasAdd_1BiasAdd(sequential_3/dense_10/MatMul_1:product:06sequential_3/dense_10/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2ђ
sequential_3/dense_10/Relu_1Relu(sequential_3/dense_10/BiasAdd_1:output:0*
T0*'
_output_shapes
:         2б
-sequential_3/dense_11/MatMul_1/ReadVariableOpReadVariableOp4sequential_3_dense_11_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0й
sequential_3/dense_11/MatMul_1MatMul*sequential_3/dense_10/Relu_1:activations:05sequential_3/dense_11/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         а
.sequential_3/dense_11/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_3_dense_11_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0Й
sequential_3/dense_11/BiasAdd_1BiasAdd(sequential_3/dense_11/MatMul_1:product:06sequential_3/dense_11/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         х
tf.stack_3/stackPack&sequential_3/dense_11/BiasAdd:output:0(sequential_3/dense_11/BiasAdd_1:output:0*
N*
T0*+
_output_shapes
:         *

axisl
IdentityIdentitytf.stack_3/stack:output:0^NoOp*
T0*+
_output_shapes
:         Ч
NoOpNoOp-^sequential_3/dense_10/BiasAdd/ReadVariableOp/^sequential_3/dense_10/BiasAdd_1/ReadVariableOp,^sequential_3/dense_10/MatMul/ReadVariableOp.^sequential_3/dense_10/MatMul_1/ReadVariableOp-^sequential_3/dense_11/BiasAdd/ReadVariableOp/^sequential_3/dense_11/BiasAdd_1/ReadVariableOp,^sequential_3/dense_11/MatMul/ReadVariableOp.^sequential_3/dense_11/MatMul_1/ReadVariableOp,^sequential_3/dense_9/BiasAdd/ReadVariableOp.^sequential_3/dense_9/BiasAdd_1/ReadVariableOp+^sequential_3/dense_9/MatMul/ReadVariableOp-^sequential_3/dense_9/MatMul_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2\
,sequential_3/dense_10/BiasAdd/ReadVariableOp,sequential_3/dense_10/BiasAdd/ReadVariableOp2`
.sequential_3/dense_10/BiasAdd_1/ReadVariableOp.sequential_3/dense_10/BiasAdd_1/ReadVariableOp2Z
+sequential_3/dense_10/MatMul/ReadVariableOp+sequential_3/dense_10/MatMul/ReadVariableOp2^
-sequential_3/dense_10/MatMul_1/ReadVariableOp-sequential_3/dense_10/MatMul_1/ReadVariableOp2\
,sequential_3/dense_11/BiasAdd/ReadVariableOp,sequential_3/dense_11/BiasAdd/ReadVariableOp2`
.sequential_3/dense_11/BiasAdd_1/ReadVariableOp.sequential_3/dense_11/BiasAdd_1/ReadVariableOp2Z
+sequential_3/dense_11/MatMul/ReadVariableOp+sequential_3/dense_11/MatMul/ReadVariableOp2^
-sequential_3/dense_11/MatMul_1/ReadVariableOp-sequential_3/dense_11/MatMul_1/ReadVariableOp2Z
+sequential_3/dense_9/BiasAdd/ReadVariableOp+sequential_3/dense_9/BiasAdd/ReadVariableOp2^
-sequential_3/dense_9/BiasAdd_1/ReadVariableOp-sequential_3/dense_9/BiasAdd_1/ReadVariableOp2X
*sequential_3/dense_9/MatMul/ReadVariableOp*sequential_3/dense_9/MatMul/ReadVariableOp2\
,sequential_3/dense_9/MatMul_1/ReadVariableOp,sequential_3/dense_9/MatMul_1/ReadVariableOp:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
К
ќ
)__inference_dense_11_layer_call_fn_173081

inputs
unknown:2
	unknown_0:
identityѕбStatefulPartitionedCallя
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8ѓ *M
fHRF
D__inference_dense_11_layer_call_and_return_conditional_losses_172409o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:         2: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         2
 
_user_specified_nameinputs
э
є
-__inference_sequential_3_layer_call_fn_172984

inputs
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityѕбStatefulPartitionedCallќ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *Q
fLRJ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172499o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
џ

З
C__inference_dense_9_layer_call_and_return_conditional_losses_173052

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:2
identityѕбBiasAdd/ReadVariableOpбMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:2*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:         2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:         2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:         : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
м
■
$__inference_signature_wrapper_172792
input_4
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityѕбStatefulPartitionedCallЗ
StatefulPartitionedCallStatefulPartitionedCallinput_4unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ **
f%R#
!__inference__wrapped_model_172358s
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*+
_output_shapes
:         `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:         
!
_user_specified_name	input_4
Ю
џ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172569
dense_9_input 
dense_9_172553:2
dense_9_172555:2!
dense_10_172558:22
dense_10_172560:2!
dense_11_172563:2
dense_11_172565:
identityѕб dense_10/StatefulPartitionedCallб dense_11/StatefulPartitionedCallбdense_9/StatefulPartitionedCallЭ
dense_9/StatefulPartitionedCallStatefulPartitionedCalldense_9_inputdense_9_172553dense_9_172555*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8ѓ *L
fGRE
C__inference_dense_9_layer_call_and_return_conditional_losses_172376Ќ
 dense_10/StatefulPartitionedCallStatefulPartitionedCall(dense_9/StatefulPartitionedCall:output:0dense_10_172558dense_10_172560*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8ѓ *M
fHRF
D__inference_dense_10_layer_call_and_return_conditional_losses_172393ў
 dense_11/StatefulPartitionedCallStatefulPartitionedCall)dense_10/StatefulPartitionedCall:output:0dense_11_172563dense_11_172565*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8ѓ *M
fHRF
D__inference_dense_11_layer_call_and_return_conditional_losses_172409x
IdentityIdentity)dense_11/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         «
NoOpNoOp!^dense_10/StatefulPartitionedCall!^dense_11/StatefulPartitionedCall ^dense_9/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2D
 dense_10/StatefulPartitionedCall dense_10/StatefulPartitionedCall2D
 dense_11/StatefulPartitionedCall dense_11/StatefulPartitionedCall2B
dense_9/StatefulPartitionedCalldense_9/StatefulPartitionedCall:V R
'
_output_shapes
:         
'
_user_specified_namedense_9_input
┼
Ћ
(__inference_dense_9_layer_call_fn_173041

inputs
unknown:2
	unknown_0:2
identityѕбStatefulPartitionedCallП
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8ѓ *L
fGRE
C__inference_dense_9_layer_call_and_return_conditional_losses_172376o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         2`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:         : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
К	
ш
D__inference_dense_11_layer_call_and_return_conditional_losses_172409

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:
identityѕбBiasAdd/ReadVariableOpбMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:2*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         _
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:         w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:         2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:         2
 
_user_specified_nameinputs
Џ

ш
D__inference_dense_10_layer_call_and_return_conditional_losses_173072

inputs0
matmul_readvariableop_resource:22-
biasadd_readvariableop_resource:2
identityѕбBiasAdd/ReadVariableOpбMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:22*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:         2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:         2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:         2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:         2
 
_user_specified_nameinputs
ф
Љ
C__inference_model_3_layer_call_and_return_conditional_losses_172738
input_4%
sequential_3_172716:2!
sequential_3_172718:2%
sequential_3_172720:22!
sequential_3_172722:2%
sequential_3_172724:2!
sequential_3_172726:
identityѕб$sequential_3/StatefulPartitionedCallб&sequential_3/StatefulPartitionedCall_1
.tf.__operators__.getitem_7/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       Ђ
0tf.__operators__.getitem_7/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        Ђ
0tf.__operators__.getitem_7/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      Л
(tf.__operators__.getitem_7/strided_sliceStridedSliceinput_47tf.__operators__.getitem_7/strided_slice/stack:output:09tf.__operators__.getitem_7/strided_slice/stack_1:output:09tf.__operators__.getitem_7/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask
.tf.__operators__.getitem_6/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        Ђ
0tf.__operators__.getitem_6/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       Ђ
0tf.__operators__.getitem_6/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      Л
(tf.__operators__.getitem_6/strided_sliceStridedSliceinput_47tf.__operators__.getitem_6/strided_slice/stack:output:09tf.__operators__.getitem_6/strided_slice/stack_1:output:09tf.__operators__.getitem_6/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_maskї
$sequential_3/StatefulPartitionedCallStatefulPartitionedCall1tf.__operators__.getitem_6/strided_slice:output:0sequential_3_172716sequential_3_172718sequential_3_172720sequential_3_172722sequential_3_172724sequential_3_172726*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *Q
fLRJ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172416ј
&sequential_3/StatefulPartitionedCall_1StatefulPartitionedCall1tf.__operators__.getitem_7/strided_slice:output:0sequential_3_172716sequential_3_172718sequential_3_172720sequential_3_172722sequential_3_172724sequential_3_172726*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *Q
fLRJ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172416├
tf.stack_3/stackPack-sequential_3/StatefulPartitionedCall:output:0/sequential_3/StatefulPartitionedCall_1:output:0*
N*
T0*+
_output_shapes
:         *

axisl
IdentityIdentitytf.stack_3/stack:output:0^NoOp*
T0*+
_output_shapes
:         ќ
NoOpNoOp%^sequential_3/StatefulPartitionedCall'^sequential_3/StatefulPartitionedCall_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2L
$sequential_3/StatefulPartitionedCall$sequential_3/StatefulPartitionedCall2P
&sequential_3/StatefulPartitionedCall_1&sequential_3/StatefulPartitionedCall_1:P L
'
_output_shapes
:         
!
_user_specified_name	input_4
ф
Љ
C__inference_model_3_layer_call_and_return_conditional_losses_172771
input_4%
sequential_3_172749:2!
sequential_3_172751:2%
sequential_3_172753:22!
sequential_3_172755:2%
sequential_3_172757:2!
sequential_3_172759:
identityѕб$sequential_3/StatefulPartitionedCallб&sequential_3/StatefulPartitionedCall_1
.tf.__operators__.getitem_7/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       Ђ
0tf.__operators__.getitem_7/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        Ђ
0tf.__operators__.getitem_7/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      Л
(tf.__operators__.getitem_7/strided_sliceStridedSliceinput_47tf.__operators__.getitem_7/strided_slice/stack:output:09tf.__operators__.getitem_7/strided_slice/stack_1:output:09tf.__operators__.getitem_7/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask
.tf.__operators__.getitem_6/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        Ђ
0tf.__operators__.getitem_6/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       Ђ
0tf.__operators__.getitem_6/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      Л
(tf.__operators__.getitem_6/strided_sliceStridedSliceinput_47tf.__operators__.getitem_6/strided_slice/stack:output:09tf.__operators__.getitem_6/strided_slice/stack_1:output:09tf.__operators__.getitem_6/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_maskї
$sequential_3/StatefulPartitionedCallStatefulPartitionedCall1tf.__operators__.getitem_6/strided_slice:output:0sequential_3_172749sequential_3_172751sequential_3_172753sequential_3_172755sequential_3_172757sequential_3_172759*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *Q
fLRJ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172499ј
&sequential_3/StatefulPartitionedCall_1StatefulPartitionedCall1tf.__operators__.getitem_7/strided_slice:output:0sequential_3_172749sequential_3_172751sequential_3_172753sequential_3_172755sequential_3_172757sequential_3_172759*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *Q
fLRJ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172499├
tf.stack_3/stackPack-sequential_3/StatefulPartitionedCall:output:0/sequential_3/StatefulPartitionedCall_1:output:0*
N*
T0*+
_output_shapes
:         *

axisl
IdentityIdentitytf.stack_3/stack:output:0^NoOp*
T0*+
_output_shapes
:         ќ
NoOpNoOp%^sequential_3/StatefulPartitionedCall'^sequential_3/StatefulPartitionedCall_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2L
$sequential_3/StatefulPartitionedCall$sequential_3/StatefulPartitionedCall2P
&sequential_3/StatefulPartitionedCall_1&sequential_3/StatefulPartitionedCall_1:P L
'
_output_shapes
:         
!
_user_specified_name	input_4
И
O
#__inference__update_step_xla_172935
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
ЮS
г	
!__inference__wrapped_model_172358
input_4M
;model_3_sequential_3_dense_9_matmul_readvariableop_resource:2J
<model_3_sequential_3_dense_9_biasadd_readvariableop_resource:2N
<model_3_sequential_3_dense_10_matmul_readvariableop_resource:22K
=model_3_sequential_3_dense_10_biasadd_readvariableop_resource:2N
<model_3_sequential_3_dense_11_matmul_readvariableop_resource:2K
=model_3_sequential_3_dense_11_biasadd_readvariableop_resource:
identityѕб4model_3/sequential_3/dense_10/BiasAdd/ReadVariableOpб6model_3/sequential_3/dense_10/BiasAdd_1/ReadVariableOpб3model_3/sequential_3/dense_10/MatMul/ReadVariableOpб5model_3/sequential_3/dense_10/MatMul_1/ReadVariableOpб4model_3/sequential_3/dense_11/BiasAdd/ReadVariableOpб6model_3/sequential_3/dense_11/BiasAdd_1/ReadVariableOpб3model_3/sequential_3/dense_11/MatMul/ReadVariableOpб5model_3/sequential_3/dense_11/MatMul_1/ReadVariableOpб3model_3/sequential_3/dense_9/BiasAdd/ReadVariableOpб5model_3/sequential_3/dense_9/BiasAdd_1/ReadVariableOpб2model_3/sequential_3/dense_9/MatMul/ReadVariableOpб4model_3/sequential_3/dense_9/MatMul_1/ReadVariableOpЄ
6model_3/tf.__operators__.getitem_7/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       Ѕ
8model_3/tf.__operators__.getitem_7/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        Ѕ
8model_3/tf.__operators__.getitem_7/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ы
0model_3/tf.__operators__.getitem_7/strided_sliceStridedSliceinput_4?model_3/tf.__operators__.getitem_7/strided_slice/stack:output:0Amodel_3/tf.__operators__.getitem_7/strided_slice/stack_1:output:0Amodel_3/tf.__operators__.getitem_7/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_maskЄ
6model_3/tf.__operators__.getitem_6/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        Ѕ
8model_3/tf.__operators__.getitem_6/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       Ѕ
8model_3/tf.__operators__.getitem_6/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ы
0model_3/tf.__operators__.getitem_6/strided_sliceStridedSliceinput_4?model_3/tf.__operators__.getitem_6/strided_slice/stack:output:0Amodel_3/tf.__operators__.getitem_6/strided_slice/stack_1:output:0Amodel_3/tf.__operators__.getitem_6/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask«
2model_3/sequential_3/dense_9/MatMul/ReadVariableOpReadVariableOp;model_3_sequential_3_dense_9_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0о
#model_3/sequential_3/dense_9/MatMulMatMul9model_3/tf.__operators__.getitem_6/strided_slice:output:0:model_3/sequential_3/dense_9/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2г
3model_3/sequential_3/dense_9/BiasAdd/ReadVariableOpReadVariableOp<model_3_sequential_3_dense_9_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0═
$model_3/sequential_3/dense_9/BiasAddBiasAdd-model_3/sequential_3/dense_9/MatMul:product:0;model_3/sequential_3/dense_9/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2і
!model_3/sequential_3/dense_9/ReluRelu-model_3/sequential_3/dense_9/BiasAdd:output:0*
T0*'
_output_shapes
:         2░
3model_3/sequential_3/dense_10/MatMul/ReadVariableOpReadVariableOp<model_3_sequential_3_dense_10_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0╬
$model_3/sequential_3/dense_10/MatMulMatMul/model_3/sequential_3/dense_9/Relu:activations:0;model_3/sequential_3/dense_10/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2«
4model_3/sequential_3/dense_10/BiasAdd/ReadVariableOpReadVariableOp=model_3_sequential_3_dense_10_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0л
%model_3/sequential_3/dense_10/BiasAddBiasAdd.model_3/sequential_3/dense_10/MatMul:product:0<model_3/sequential_3/dense_10/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2ї
"model_3/sequential_3/dense_10/ReluRelu.model_3/sequential_3/dense_10/BiasAdd:output:0*
T0*'
_output_shapes
:         2░
3model_3/sequential_3/dense_11/MatMul/ReadVariableOpReadVariableOp<model_3_sequential_3_dense_11_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0¤
$model_3/sequential_3/dense_11/MatMulMatMul0model_3/sequential_3/dense_10/Relu:activations:0;model_3/sequential_3/dense_11/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         «
4model_3/sequential_3/dense_11/BiasAdd/ReadVariableOpReadVariableOp=model_3_sequential_3_dense_11_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0л
%model_3/sequential_3/dense_11/BiasAddBiasAdd.model_3/sequential_3/dense_11/MatMul:product:0<model_3/sequential_3/dense_11/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         ░
4model_3/sequential_3/dense_9/MatMul_1/ReadVariableOpReadVariableOp;model_3_sequential_3_dense_9_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0┌
%model_3/sequential_3/dense_9/MatMul_1MatMul9model_3/tf.__operators__.getitem_7/strided_slice:output:0<model_3/sequential_3/dense_9/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2«
5model_3/sequential_3/dense_9/BiasAdd_1/ReadVariableOpReadVariableOp<model_3_sequential_3_dense_9_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0М
&model_3/sequential_3/dense_9/BiasAdd_1BiasAdd/model_3/sequential_3/dense_9/MatMul_1:product:0=model_3/sequential_3/dense_9/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2ј
#model_3/sequential_3/dense_9/Relu_1Relu/model_3/sequential_3/dense_9/BiasAdd_1:output:0*
T0*'
_output_shapes
:         2▓
5model_3/sequential_3/dense_10/MatMul_1/ReadVariableOpReadVariableOp<model_3_sequential_3_dense_10_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0н
&model_3/sequential_3/dense_10/MatMul_1MatMul1model_3/sequential_3/dense_9/Relu_1:activations:0=model_3/sequential_3/dense_10/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2░
6model_3/sequential_3/dense_10/BiasAdd_1/ReadVariableOpReadVariableOp=model_3_sequential_3_dense_10_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0о
'model_3/sequential_3/dense_10/BiasAdd_1BiasAdd0model_3/sequential_3/dense_10/MatMul_1:product:0>model_3/sequential_3/dense_10/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2љ
$model_3/sequential_3/dense_10/Relu_1Relu0model_3/sequential_3/dense_10/BiasAdd_1:output:0*
T0*'
_output_shapes
:         2▓
5model_3/sequential_3/dense_11/MatMul_1/ReadVariableOpReadVariableOp<model_3_sequential_3_dense_11_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0Н
&model_3/sequential_3/dense_11/MatMul_1MatMul2model_3/sequential_3/dense_10/Relu_1:activations:0=model_3/sequential_3/dense_11/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         ░
6model_3/sequential_3/dense_11/BiasAdd_1/ReadVariableOpReadVariableOp=model_3_sequential_3_dense_11_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0о
'model_3/sequential_3/dense_11/BiasAdd_1BiasAdd0model_3/sequential_3/dense_11/MatMul_1:product:0>model_3/sequential_3/dense_11/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         ═
model_3/tf.stack_3/stackPack.model_3/sequential_3/dense_11/BiasAdd:output:00model_3/sequential_3/dense_11/BiasAdd_1:output:0*
N*
T0*+
_output_shapes
:         *

axist
IdentityIdentity!model_3/tf.stack_3/stack:output:0^NoOp*
T0*+
_output_shapes
:         ▄
NoOpNoOp5^model_3/sequential_3/dense_10/BiasAdd/ReadVariableOp7^model_3/sequential_3/dense_10/BiasAdd_1/ReadVariableOp4^model_3/sequential_3/dense_10/MatMul/ReadVariableOp6^model_3/sequential_3/dense_10/MatMul_1/ReadVariableOp5^model_3/sequential_3/dense_11/BiasAdd/ReadVariableOp7^model_3/sequential_3/dense_11/BiasAdd_1/ReadVariableOp4^model_3/sequential_3/dense_11/MatMul/ReadVariableOp6^model_3/sequential_3/dense_11/MatMul_1/ReadVariableOp4^model_3/sequential_3/dense_9/BiasAdd/ReadVariableOp6^model_3/sequential_3/dense_9/BiasAdd_1/ReadVariableOp3^model_3/sequential_3/dense_9/MatMul/ReadVariableOp5^model_3/sequential_3/dense_9/MatMul_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2l
4model_3/sequential_3/dense_10/BiasAdd/ReadVariableOp4model_3/sequential_3/dense_10/BiasAdd/ReadVariableOp2p
6model_3/sequential_3/dense_10/BiasAdd_1/ReadVariableOp6model_3/sequential_3/dense_10/BiasAdd_1/ReadVariableOp2j
3model_3/sequential_3/dense_10/MatMul/ReadVariableOp3model_3/sequential_3/dense_10/MatMul/ReadVariableOp2n
5model_3/sequential_3/dense_10/MatMul_1/ReadVariableOp5model_3/sequential_3/dense_10/MatMul_1/ReadVariableOp2l
4model_3/sequential_3/dense_11/BiasAdd/ReadVariableOp4model_3/sequential_3/dense_11/BiasAdd/ReadVariableOp2p
6model_3/sequential_3/dense_11/BiasAdd_1/ReadVariableOp6model_3/sequential_3/dense_11/BiasAdd_1/ReadVariableOp2j
3model_3/sequential_3/dense_11/MatMul/ReadVariableOp3model_3/sequential_3/dense_11/MatMul/ReadVariableOp2n
5model_3/sequential_3/dense_11/MatMul_1/ReadVariableOp5model_3/sequential_3/dense_11/MatMul_1/ReadVariableOp2j
3model_3/sequential_3/dense_9/BiasAdd/ReadVariableOp3model_3/sequential_3/dense_9/BiasAdd/ReadVariableOp2n
5model_3/sequential_3/dense_9/BiasAdd_1/ReadVariableOp5model_3/sequential_3/dense_9/BiasAdd_1/ReadVariableOp2h
2model_3/sequential_3/dense_9/MatMul/ReadVariableOp2model_3/sequential_3/dense_9/MatMul/ReadVariableOp2l
4model_3/sequential_3/dense_9/MatMul_1/ReadVariableOp4model_3/sequential_3/dense_9/MatMul_1/ReadVariableOp:P L
'
_output_shapes
:         
!
_user_specified_name	input_4
Џ

ш
D__inference_dense_10_layer_call_and_return_conditional_losses_172393

inputs0
matmul_readvariableop_resource:22-
biasadd_readvariableop_resource:2
identityѕбBiasAdd/ReadVariableOpбMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:22*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:         2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:         2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:         2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:         2
 
_user_specified_nameinputs
ш
Ђ
(__inference_model_3_layer_call_fn_172826

inputs
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityѕбStatefulPartitionedCallЋ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *L
fGRE
C__inference_model_3_layer_call_and_return_conditional_losses_172673s
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*+
_output_shapes
:         `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
г
K
#__inference__update_step_xla_172940
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
г
K
#__inference__update_step_xla_172930
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
ѕ
Њ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172416

inputs 
dense_9_172377:2
dense_9_172379:2!
dense_10_172394:22
dense_10_172396:2!
dense_11_172410:2
dense_11_172412:
identityѕб dense_10/StatefulPartitionedCallб dense_11/StatefulPartitionedCallбdense_9/StatefulPartitionedCallы
dense_9/StatefulPartitionedCallStatefulPartitionedCallinputsdense_9_172377dense_9_172379*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8ѓ *L
fGRE
C__inference_dense_9_layer_call_and_return_conditional_losses_172376Ќ
 dense_10/StatefulPartitionedCallStatefulPartitionedCall(dense_9/StatefulPartitionedCall:output:0dense_10_172394dense_10_172396*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8ѓ *M
fHRF
D__inference_dense_10_layer_call_and_return_conditional_losses_172393ў
 dense_11/StatefulPartitionedCallStatefulPartitionedCall)dense_10/StatefulPartitionedCall:output:0dense_11_172410dense_11_172412*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8ѓ *M
fHRF
D__inference_dense_11_layer_call_and_return_conditional_losses_172409x
IdentityIdentity)dense_11/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         «
NoOpNoOp!^dense_10/StatefulPartitionedCall!^dense_11/StatefulPartitionedCall ^dense_9/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2D
 dense_10/StatefulPartitionedCall dense_10/StatefulPartitionedCall2D
 dense_11/StatefulPartitionedCall dense_11/StatefulPartitionedCall2B
dense_9/StatefulPartitionedCalldense_9/StatefulPartitionedCall:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
И
O
#__inference__update_step_xla_172925
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
оJ
й
C__inference_model_3_layer_call_and_return_conditional_losses_172920

inputsE
3sequential_3_dense_9_matmul_readvariableop_resource:2B
4sequential_3_dense_9_biasadd_readvariableop_resource:2F
4sequential_3_dense_10_matmul_readvariableop_resource:22C
5sequential_3_dense_10_biasadd_readvariableop_resource:2F
4sequential_3_dense_11_matmul_readvariableop_resource:2C
5sequential_3_dense_11_biasadd_readvariableop_resource:
identityѕб,sequential_3/dense_10/BiasAdd/ReadVariableOpб.sequential_3/dense_10/BiasAdd_1/ReadVariableOpб+sequential_3/dense_10/MatMul/ReadVariableOpб-sequential_3/dense_10/MatMul_1/ReadVariableOpб,sequential_3/dense_11/BiasAdd/ReadVariableOpб.sequential_3/dense_11/BiasAdd_1/ReadVariableOpб+sequential_3/dense_11/MatMul/ReadVariableOpб-sequential_3/dense_11/MatMul_1/ReadVariableOpб+sequential_3/dense_9/BiasAdd/ReadVariableOpб-sequential_3/dense_9/BiasAdd_1/ReadVariableOpб*sequential_3/dense_9/MatMul/ReadVariableOpб,sequential_3/dense_9/MatMul_1/ReadVariableOp
.tf.__operators__.getitem_7/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       Ђ
0tf.__operators__.getitem_7/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        Ђ
0tf.__operators__.getitem_7/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      л
(tf.__operators__.getitem_7/strided_sliceStridedSliceinputs7tf.__operators__.getitem_7/strided_slice/stack:output:09tf.__operators__.getitem_7/strided_slice/stack_1:output:09tf.__operators__.getitem_7/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask
.tf.__operators__.getitem_6/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        Ђ
0tf.__operators__.getitem_6/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       Ђ
0tf.__operators__.getitem_6/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      л
(tf.__operators__.getitem_6/strided_sliceStridedSliceinputs7tf.__operators__.getitem_6/strided_slice/stack:output:09tf.__operators__.getitem_6/strided_slice/stack_1:output:09tf.__operators__.getitem_6/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_maskъ
*sequential_3/dense_9/MatMul/ReadVariableOpReadVariableOp3sequential_3_dense_9_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0Й
sequential_3/dense_9/MatMulMatMul1tf.__operators__.getitem_6/strided_slice:output:02sequential_3/dense_9/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2ю
+sequential_3/dense_9/BiasAdd/ReadVariableOpReadVariableOp4sequential_3_dense_9_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0х
sequential_3/dense_9/BiasAddBiasAdd%sequential_3/dense_9/MatMul:product:03sequential_3/dense_9/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2z
sequential_3/dense_9/ReluRelu%sequential_3/dense_9/BiasAdd:output:0*
T0*'
_output_shapes
:         2а
+sequential_3/dense_10/MatMul/ReadVariableOpReadVariableOp4sequential_3_dense_10_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Х
sequential_3/dense_10/MatMulMatMul'sequential_3/dense_9/Relu:activations:03sequential_3/dense_10/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2ъ
,sequential_3/dense_10/BiasAdd/ReadVariableOpReadVariableOp5sequential_3_dense_10_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0И
sequential_3/dense_10/BiasAddBiasAdd&sequential_3/dense_10/MatMul:product:04sequential_3/dense_10/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2|
sequential_3/dense_10/ReluRelu&sequential_3/dense_10/BiasAdd:output:0*
T0*'
_output_shapes
:         2а
+sequential_3/dense_11/MatMul/ReadVariableOpReadVariableOp4sequential_3_dense_11_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0и
sequential_3/dense_11/MatMulMatMul(sequential_3/dense_10/Relu:activations:03sequential_3/dense_11/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         ъ
,sequential_3/dense_11/BiasAdd/ReadVariableOpReadVariableOp5sequential_3_dense_11_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0И
sequential_3/dense_11/BiasAddBiasAdd&sequential_3/dense_11/MatMul:product:04sequential_3/dense_11/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         а
,sequential_3/dense_9/MatMul_1/ReadVariableOpReadVariableOp3sequential_3_dense_9_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0┬
sequential_3/dense_9/MatMul_1MatMul1tf.__operators__.getitem_7/strided_slice:output:04sequential_3/dense_9/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2ъ
-sequential_3/dense_9/BiasAdd_1/ReadVariableOpReadVariableOp4sequential_3_dense_9_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0╗
sequential_3/dense_9/BiasAdd_1BiasAdd'sequential_3/dense_9/MatMul_1:product:05sequential_3/dense_9/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2~
sequential_3/dense_9/Relu_1Relu'sequential_3/dense_9/BiasAdd_1:output:0*
T0*'
_output_shapes
:         2б
-sequential_3/dense_10/MatMul_1/ReadVariableOpReadVariableOp4sequential_3_dense_10_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0╝
sequential_3/dense_10/MatMul_1MatMul)sequential_3/dense_9/Relu_1:activations:05sequential_3/dense_10/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2а
.sequential_3/dense_10/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_3_dense_10_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0Й
sequential_3/dense_10/BiasAdd_1BiasAdd(sequential_3/dense_10/MatMul_1:product:06sequential_3/dense_10/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2ђ
sequential_3/dense_10/Relu_1Relu(sequential_3/dense_10/BiasAdd_1:output:0*
T0*'
_output_shapes
:         2б
-sequential_3/dense_11/MatMul_1/ReadVariableOpReadVariableOp4sequential_3_dense_11_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0й
sequential_3/dense_11/MatMul_1MatMul*sequential_3/dense_10/Relu_1:activations:05sequential_3/dense_11/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         а
.sequential_3/dense_11/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_3_dense_11_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0Й
sequential_3/dense_11/BiasAdd_1BiasAdd(sequential_3/dense_11/MatMul_1:product:06sequential_3/dense_11/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         х
tf.stack_3/stackPack&sequential_3/dense_11/BiasAdd:output:0(sequential_3/dense_11/BiasAdd_1:output:0*
N*
T0*+
_output_shapes
:         *

axisl
IdentityIdentitytf.stack_3/stack:output:0^NoOp*
T0*+
_output_shapes
:         Ч
NoOpNoOp-^sequential_3/dense_10/BiasAdd/ReadVariableOp/^sequential_3/dense_10/BiasAdd_1/ReadVariableOp,^sequential_3/dense_10/MatMul/ReadVariableOp.^sequential_3/dense_10/MatMul_1/ReadVariableOp-^sequential_3/dense_11/BiasAdd/ReadVariableOp/^sequential_3/dense_11/BiasAdd_1/ReadVariableOp,^sequential_3/dense_11/MatMul/ReadVariableOp.^sequential_3/dense_11/MatMul_1/ReadVariableOp,^sequential_3/dense_9/BiasAdd/ReadVariableOp.^sequential_3/dense_9/BiasAdd_1/ReadVariableOp+^sequential_3/dense_9/MatMul/ReadVariableOp-^sequential_3/dense_9/MatMul_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2\
,sequential_3/dense_10/BiasAdd/ReadVariableOp,sequential_3/dense_10/BiasAdd/ReadVariableOp2`
.sequential_3/dense_10/BiasAdd_1/ReadVariableOp.sequential_3/dense_10/BiasAdd_1/ReadVariableOp2Z
+sequential_3/dense_10/MatMul/ReadVariableOp+sequential_3/dense_10/MatMul/ReadVariableOp2^
-sequential_3/dense_10/MatMul_1/ReadVariableOp-sequential_3/dense_10/MatMul_1/ReadVariableOp2\
,sequential_3/dense_11/BiasAdd/ReadVariableOp,sequential_3/dense_11/BiasAdd/ReadVariableOp2`
.sequential_3/dense_11/BiasAdd_1/ReadVariableOp.sequential_3/dense_11/BiasAdd_1/ReadVariableOp2Z
+sequential_3/dense_11/MatMul/ReadVariableOp+sequential_3/dense_11/MatMul/ReadVariableOp2^
-sequential_3/dense_11/MatMul_1/ReadVariableOp-sequential_3/dense_11/MatMul_1/ReadVariableOp2Z
+sequential_3/dense_9/BiasAdd/ReadVariableOp+sequential_3/dense_9/BiasAdd/ReadVariableOp2^
-sequential_3/dense_9/BiasAdd_1/ReadVariableOp-sequential_3/dense_9/BiasAdd_1/ReadVariableOp2X
*sequential_3/dense_9/MatMul/ReadVariableOp*sequential_3/dense_9/MatMul/ReadVariableOp2\
,sequential_3/dense_9/MatMul_1/ReadVariableOp,sequential_3/dense_9/MatMul_1/ReadVariableOp:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
ш
Ђ
(__inference_model_3_layer_call_fn_172809

inputs
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityѕбStatefulPartitionedCallЋ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *L
fGRE
C__inference_model_3_layer_call_and_return_conditional_losses_172606s
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*+
_output_shapes
:         `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
ї	
Ї
-__inference_sequential_3_layer_call_fn_172431
dense_9_input
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityѕбStatefulPartitionedCallЮ
StatefulPartitionedCallStatefulPartitionedCalldense_9_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *Q
fLRJ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172416o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:V R
'
_output_shapes
:         
'
_user_specified_namedense_9_input
ї	
Ї
-__inference_sequential_3_layer_call_fn_172531
dense_9_input
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityѕбStatefulPartitionedCallЮ
StatefulPartitionedCallStatefulPartitionedCalldense_9_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *Q
fLRJ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172499o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:V R
'
_output_shapes
:         
'
_user_specified_namedense_9_input
Э
ѓ
(__inference_model_3_layer_call_fn_172621
input_4
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityѕбStatefulPartitionedCallќ
StatefulPartitionedCallStatefulPartitionedCallinput_4unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *L
fGRE
C__inference_model_3_layer_call_and_return_conditional_losses_172606s
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*+
_output_shapes
:         `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:         
!
_user_specified_name	input_4
э
є
-__inference_sequential_3_layer_call_fn_172967

inputs
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityѕбStatefulPartitionedCallќ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *Q
fLRJ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172416o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
д
љ
C__inference_model_3_layer_call_and_return_conditional_losses_172673

inputs%
sequential_3_172651:2!
sequential_3_172653:2%
sequential_3_172655:22!
sequential_3_172657:2%
sequential_3_172659:2!
sequential_3_172661:
identityѕб$sequential_3/StatefulPartitionedCallб&sequential_3/StatefulPartitionedCall_1
.tf.__operators__.getitem_7/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       Ђ
0tf.__operators__.getitem_7/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        Ђ
0tf.__operators__.getitem_7/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      л
(tf.__operators__.getitem_7/strided_sliceStridedSliceinputs7tf.__operators__.getitem_7/strided_slice/stack:output:09tf.__operators__.getitem_7/strided_slice/stack_1:output:09tf.__operators__.getitem_7/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask
.tf.__operators__.getitem_6/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        Ђ
0tf.__operators__.getitem_6/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       Ђ
0tf.__operators__.getitem_6/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      л
(tf.__operators__.getitem_6/strided_sliceStridedSliceinputs7tf.__operators__.getitem_6/strided_slice/stack:output:09tf.__operators__.getitem_6/strided_slice/stack_1:output:09tf.__operators__.getitem_6/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_maskї
$sequential_3/StatefulPartitionedCallStatefulPartitionedCall1tf.__operators__.getitem_6/strided_slice:output:0sequential_3_172651sequential_3_172653sequential_3_172655sequential_3_172657sequential_3_172659sequential_3_172661*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *Q
fLRJ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172499ј
&sequential_3/StatefulPartitionedCall_1StatefulPartitionedCall1tf.__operators__.getitem_7/strided_slice:output:0sequential_3_172651sequential_3_172653sequential_3_172655sequential_3_172657sequential_3_172659sequential_3_172661*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *Q
fLRJ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172499├
tf.stack_3/stackPack-sequential_3/StatefulPartitionedCall:output:0/sequential_3/StatefulPartitionedCall_1:output:0*
N*
T0*+
_output_shapes
:         *

axisl
IdentityIdentitytf.stack_3/stack:output:0^NoOp*
T0*+
_output_shapes
:         ќ
NoOpNoOp%^sequential_3/StatefulPartitionedCall'^sequential_3/StatefulPartitionedCall_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2L
$sequential_3/StatefulPartitionedCall$sequential_3/StatefulPartitionedCall2P
&sequential_3/StatefulPartitionedCall_1&sequential_3/StatefulPartitionedCall_1:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
▄4
Ч	
__inference__traced_save_173186
file_prefix-
)savev2_dense_9_kernel_read_readvariableop+
'savev2_dense_9_bias_read_readvariableop.
*savev2_dense_10_kernel_read_readvariableop,
(savev2_dense_10_bias_read_readvariableop.
*savev2_dense_11_kernel_read_readvariableop,
(savev2_dense_11_bias_read_readvariableop(
$savev2_iteration_read_readvariableop	,
(savev2_learning_rate_read_readvariableop4
0savev2_adam_m_dense_9_kernel_read_readvariableop4
0savev2_adam_v_dense_9_kernel_read_readvariableop2
.savev2_adam_m_dense_9_bias_read_readvariableop2
.savev2_adam_v_dense_9_bias_read_readvariableop5
1savev2_adam_m_dense_10_kernel_read_readvariableop5
1savev2_adam_v_dense_10_kernel_read_readvariableop3
/savev2_adam_m_dense_10_bias_read_readvariableop3
/savev2_adam_v_dense_10_bias_read_readvariableop5
1savev2_adam_m_dense_11_kernel_read_readvariableop5
1savev2_adam_v_dense_11_kernel_read_readvariableop3
/savev2_adam_m_dense_11_bias_read_readvariableop3
/savev2_adam_v_dense_11_bias_read_readvariableop&
"savev2_total_1_read_readvariableop&
"savev2_count_1_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop
savev2_const

identity_1ѕбMergeV2Checkpointsw
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
_temp/partЂ
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
value	B : Њ
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: а

SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*╔	
value┐	B╝	B&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUEB3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPHЪ
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*E
value<B:B B B B B B B B B B B B B B B B B B B B B B B B B а

SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0)savev2_dense_9_kernel_read_readvariableop'savev2_dense_9_bias_read_readvariableop*savev2_dense_10_kernel_read_readvariableop(savev2_dense_10_bias_read_readvariableop*savev2_dense_11_kernel_read_readvariableop(savev2_dense_11_bias_read_readvariableop$savev2_iteration_read_readvariableop(savev2_learning_rate_read_readvariableop0savev2_adam_m_dense_9_kernel_read_readvariableop0savev2_adam_v_dense_9_kernel_read_readvariableop.savev2_adam_m_dense_9_bias_read_readvariableop.savev2_adam_v_dense_9_bias_read_readvariableop1savev2_adam_m_dense_10_kernel_read_readvariableop1savev2_adam_v_dense_10_kernel_read_readvariableop/savev2_adam_m_dense_10_bias_read_readvariableop/savev2_adam_v_dense_10_bias_read_readvariableop1savev2_adam_m_dense_11_kernel_read_readvariableop1savev2_adam_v_dense_11_kernel_read_readvariableop/savev2_adam_m_dense_11_bias_read_readvariableop/savev2_adam_v_dense_11_bias_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableopsavev2_const"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *'
dtypes
2	љ
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:│
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

identity_1Identity_1:output:0*х
_input_shapesБ
а: :2:2:22:2:2:: : :2:2:2:2:22:22:2:2:2:2::: : : : : 2(
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
Ю
џ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172550
dense_9_input 
dense_9_172534:2
dense_9_172536:2!
dense_10_172539:22
dense_10_172541:2!
dense_11_172544:2
dense_11_172546:
identityѕб dense_10/StatefulPartitionedCallб dense_11/StatefulPartitionedCallбdense_9/StatefulPartitionedCallЭ
dense_9/StatefulPartitionedCallStatefulPartitionedCalldense_9_inputdense_9_172534dense_9_172536*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8ѓ *L
fGRE
C__inference_dense_9_layer_call_and_return_conditional_losses_172376Ќ
 dense_10/StatefulPartitionedCallStatefulPartitionedCall(dense_9/StatefulPartitionedCall:output:0dense_10_172539dense_10_172541*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8ѓ *M
fHRF
D__inference_dense_10_layer_call_and_return_conditional_losses_172393ў
 dense_11/StatefulPartitionedCallStatefulPartitionedCall)dense_10/StatefulPartitionedCall:output:0dense_11_172544dense_11_172546*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8ѓ *M
fHRF
D__inference_dense_11_layer_call_and_return_conditional_losses_172409x
IdentityIdentity)dense_11/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         «
NoOpNoOp!^dense_10/StatefulPartitionedCall!^dense_11/StatefulPartitionedCall ^dense_9/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2D
 dense_10/StatefulPartitionedCall dense_10/StatefulPartitionedCall2D
 dense_11/StatefulPartitionedCall dense_11/StatefulPartitionedCall2B
dense_9/StatefulPartitionedCalldense_9/StatefulPartitionedCall:V R
'
_output_shapes
:         
'
_user_specified_namedense_9_input
Е
Ё
H__inference_sequential_3_layer_call_and_return_conditional_losses_173008

inputs8
&dense_9_matmul_readvariableop_resource:25
'dense_9_biasadd_readvariableop_resource:29
'dense_10_matmul_readvariableop_resource:226
(dense_10_biasadd_readvariableop_resource:29
'dense_11_matmul_readvariableop_resource:26
(dense_11_biasadd_readvariableop_resource:
identityѕбdense_10/BiasAdd/ReadVariableOpбdense_10/MatMul/ReadVariableOpбdense_11/BiasAdd/ReadVariableOpбdense_11/MatMul/ReadVariableOpбdense_9/BiasAdd/ReadVariableOpбdense_9/MatMul/ReadVariableOpё
dense_9/MatMul/ReadVariableOpReadVariableOp&dense_9_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0y
dense_9/MatMulMatMulinputs%dense_9/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2ѓ
dense_9/BiasAdd/ReadVariableOpReadVariableOp'dense_9_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0ј
dense_9/BiasAddBiasAdddense_9/MatMul:product:0&dense_9/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2`
dense_9/ReluReludense_9/BiasAdd:output:0*
T0*'
_output_shapes
:         2є
dense_10/MatMul/ReadVariableOpReadVariableOp'dense_10_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Ј
dense_10/MatMulMatMuldense_9/Relu:activations:0&dense_10/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2ё
dense_10/BiasAdd/ReadVariableOpReadVariableOp(dense_10_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0Љ
dense_10/BiasAddBiasAdddense_10/MatMul:product:0'dense_10/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2b
dense_10/ReluReludense_10/BiasAdd:output:0*
T0*'
_output_shapes
:         2є
dense_11/MatMul/ReadVariableOpReadVariableOp'dense_11_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0љ
dense_11/MatMulMatMuldense_10/Relu:activations:0&dense_11/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         ё
dense_11/BiasAdd/ReadVariableOpReadVariableOp(dense_11_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0Љ
dense_11/BiasAddBiasAdddense_11/MatMul:product:0'dense_11/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         h
IdentityIdentitydense_11/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:         Ї
NoOpNoOp ^dense_10/BiasAdd/ReadVariableOp^dense_10/MatMul/ReadVariableOp ^dense_11/BiasAdd/ReadVariableOp^dense_11/MatMul/ReadVariableOp^dense_9/BiasAdd/ReadVariableOp^dense_9/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2B
dense_10/BiasAdd/ReadVariableOpdense_10/BiasAdd/ReadVariableOp2@
dense_10/MatMul/ReadVariableOpdense_10/MatMul/ReadVariableOp2B
dense_11/BiasAdd/ReadVariableOpdense_11/BiasAdd/ReadVariableOp2@
dense_11/MatMul/ReadVariableOpdense_11/MatMul/ReadVariableOp2@
dense_9/BiasAdd/ReadVariableOpdense_9/BiasAdd/ReadVariableOp2>
dense_9/MatMul/ReadVariableOpdense_9/MatMul/ReadVariableOp:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
д
љ
C__inference_model_3_layer_call_and_return_conditional_losses_172606

inputs%
sequential_3_172584:2!
sequential_3_172586:2%
sequential_3_172588:22!
sequential_3_172590:2%
sequential_3_172592:2!
sequential_3_172594:
identityѕб$sequential_3/StatefulPartitionedCallб&sequential_3/StatefulPartitionedCall_1
.tf.__operators__.getitem_7/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"       Ђ
0tf.__operators__.getitem_7/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        Ђ
0tf.__operators__.getitem_7/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      л
(tf.__operators__.getitem_7/strided_sliceStridedSliceinputs7tf.__operators__.getitem_7/strided_slice/stack:output:09tf.__operators__.getitem_7/strided_slice/stack_1:output:09tf.__operators__.getitem_7/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask
.tf.__operators__.getitem_6/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        Ђ
0tf.__operators__.getitem_6/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       Ђ
0tf.__operators__.getitem_6/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      л
(tf.__operators__.getitem_6/strided_sliceStridedSliceinputs7tf.__operators__.getitem_6/strided_slice/stack:output:09tf.__operators__.getitem_6/strided_slice/stack_1:output:09tf.__operators__.getitem_6/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_maskї
$sequential_3/StatefulPartitionedCallStatefulPartitionedCall1tf.__operators__.getitem_6/strided_slice:output:0sequential_3_172584sequential_3_172586sequential_3_172588sequential_3_172590sequential_3_172592sequential_3_172594*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *Q
fLRJ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172416ј
&sequential_3/StatefulPartitionedCall_1StatefulPartitionedCall1tf.__operators__.getitem_7/strided_slice:output:0sequential_3_172584sequential_3_172586sequential_3_172588sequential_3_172590sequential_3_172592sequential_3_172594*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *Q
fLRJ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172416├
tf.stack_3/stackPack-sequential_3/StatefulPartitionedCall:output:0/sequential_3/StatefulPartitionedCall_1:output:0*
N*
T0*+
_output_shapes
:         *

axisl
IdentityIdentitytf.stack_3/stack:output:0^NoOp*
T0*+
_output_shapes
:         ќ
NoOpNoOp%^sequential_3/StatefulPartitionedCall'^sequential_3/StatefulPartitionedCall_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2L
$sequential_3/StatefulPartitionedCall$sequential_3/StatefulPartitionedCall2P
&sequential_3/StatefulPartitionedCall_1&sequential_3/StatefulPartitionedCall_1:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
К	
ш
D__inference_dense_11_layer_call_and_return_conditional_losses_173091

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:
identityѕбBiasAdd/ReadVariableOpбMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:2*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         _
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:         w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:         2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:         2
 
_user_specified_nameinputs
ѕ
Њ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172499

inputs 
dense_9_172483:2
dense_9_172485:2!
dense_10_172488:22
dense_10_172490:2!
dense_11_172493:2
dense_11_172495:
identityѕб dense_10/StatefulPartitionedCallб dense_11/StatefulPartitionedCallбdense_9/StatefulPartitionedCallы
dense_9/StatefulPartitionedCallStatefulPartitionedCallinputsdense_9_172483dense_9_172485*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8ѓ *L
fGRE
C__inference_dense_9_layer_call_and_return_conditional_losses_172376Ќ
 dense_10/StatefulPartitionedCallStatefulPartitionedCall(dense_9/StatefulPartitionedCall:output:0dense_10_172488dense_10_172490*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8ѓ *M
fHRF
D__inference_dense_10_layer_call_and_return_conditional_losses_172393ў
 dense_11/StatefulPartitionedCallStatefulPartitionedCall)dense_10/StatefulPartitionedCall:output:0dense_11_172493dense_11_172495*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8ѓ *M
fHRF
D__inference_dense_11_layer_call_and_return_conditional_losses_172409x
IdentityIdentity)dense_11/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         «
NoOpNoOp!^dense_10/StatefulPartitionedCall!^dense_11/StatefulPartitionedCall ^dense_9/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2D
 dense_10/StatefulPartitionedCall dense_10/StatefulPartitionedCall2D
 dense_11/StatefulPartitionedCall dense_11/StatefulPartitionedCall2B
dense_9/StatefulPartitionedCalldense_9/StatefulPartitionedCall:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
Э
ѓ
(__inference_model_3_layer_call_fn_172705
input_4
unknown:2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityѕбStatefulPartitionedCallќ
StatefulPartitionedCallStatefulPartitionedCallinput_4unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *+
_output_shapes
:         *(
_read_only_resource_inputs

*2
config_proto" 

CPU

GPU2 *0J 8ѓ *L
fGRE
C__inference_model_3_layer_call_and_return_conditional_losses_172673s
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*+
_output_shapes
:         `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:         
!
_user_specified_name	input_4
Ёf
Ц
"__inference__traced_restore_173268
file_prefix1
assignvariableop_dense_9_kernel:2-
assignvariableop_1_dense_9_bias:24
"assignvariableop_2_dense_10_kernel:22.
 assignvariableop_3_dense_10_bias:24
"assignvariableop_4_dense_11_kernel:2.
 assignvariableop_5_dense_11_bias:&
assignvariableop_6_iteration:	 *
 assignvariableop_7_learning_rate: :
(assignvariableop_8_adam_m_dense_9_kernel:2:
(assignvariableop_9_adam_v_dense_9_kernel:25
'assignvariableop_10_adam_m_dense_9_bias:25
'assignvariableop_11_adam_v_dense_9_bias:2<
*assignvariableop_12_adam_m_dense_10_kernel:22<
*assignvariableop_13_adam_v_dense_10_kernel:226
(assignvariableop_14_adam_m_dense_10_bias:26
(assignvariableop_15_adam_v_dense_10_bias:2<
*assignvariableop_16_adam_m_dense_11_kernel:2<
*assignvariableop_17_adam_v_dense_11_kernel:26
(assignvariableop_18_adam_m_dense_11_bias:6
(assignvariableop_19_adam_v_dense_11_bias:%
assignvariableop_20_total_1: %
assignvariableop_21_count_1: #
assignvariableop_22_total: #
assignvariableop_23_count: 
identity_25ѕбAssignVariableOpбAssignVariableOp_1бAssignVariableOp_10бAssignVariableOp_11бAssignVariableOp_12бAssignVariableOp_13бAssignVariableOp_14бAssignVariableOp_15бAssignVariableOp_16бAssignVariableOp_17бAssignVariableOp_18бAssignVariableOp_19бAssignVariableOp_2бAssignVariableOp_20бAssignVariableOp_21бAssignVariableOp_22бAssignVariableOp_23бAssignVariableOp_3бAssignVariableOp_4бAssignVariableOp_5бAssignVariableOp_6бAssignVariableOp_7бAssignVariableOp_8бAssignVariableOp_9Б

RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*╔	
value┐	B╝	B&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUEB3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPHб
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*E
value<B:B B B B B B B B B B B B B B B B B B B B B B B B B Џ
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*x
_output_shapesf
d:::::::::::::::::::::::::*'
dtypes
2	[
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:▓
AssignVariableOpAssignVariableOpassignvariableop_dense_9_kernelIdentity:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:Х
AssignVariableOp_1AssignVariableOpassignvariableop_1_dense_9_biasIdentity_1:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:╣
AssignVariableOp_2AssignVariableOp"assignvariableop_2_dense_10_kernelIdentity_2:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:и
AssignVariableOp_3AssignVariableOp assignvariableop_3_dense_10_biasIdentity_3:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:╣
AssignVariableOp_4AssignVariableOp"assignvariableop_4_dense_11_kernelIdentity_4:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:и
AssignVariableOp_5AssignVariableOp assignvariableop_5_dense_11_biasIdentity_5:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0	*
_output_shapes
:│
AssignVariableOp_6AssignVariableOpassignvariableop_6_iterationIdentity_6:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0	]

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:и
AssignVariableOp_7AssignVariableOp assignvariableop_7_learning_rateIdentity_7:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:┐
AssignVariableOp_8AssignVariableOp(assignvariableop_8_adam_m_dense_9_kernelIdentity_8:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:┐
AssignVariableOp_9AssignVariableOp(assignvariableop_9_adam_v_dense_9_kernelIdentity_9:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:└
AssignVariableOp_10AssignVariableOp'assignvariableop_10_adam_m_dense_9_biasIdentity_10:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:└
AssignVariableOp_11AssignVariableOp'assignvariableop_11_adam_v_dense_9_biasIdentity_11:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:├
AssignVariableOp_12AssignVariableOp*assignvariableop_12_adam_m_dense_10_kernelIdentity_12:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:├
AssignVariableOp_13AssignVariableOp*assignvariableop_13_adam_v_dense_10_kernelIdentity_13:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:┴
AssignVariableOp_14AssignVariableOp(assignvariableop_14_adam_m_dense_10_biasIdentity_14:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:┴
AssignVariableOp_15AssignVariableOp(assignvariableop_15_adam_v_dense_10_biasIdentity_15:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:├
AssignVariableOp_16AssignVariableOp*assignvariableop_16_adam_m_dense_11_kernelIdentity_16:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:├
AssignVariableOp_17AssignVariableOp*assignvariableop_17_adam_v_dense_11_kernelIdentity_17:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:┴
AssignVariableOp_18AssignVariableOp(assignvariableop_18_adam_m_dense_11_biasIdentity_18:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:┴
AssignVariableOp_19AssignVariableOp(assignvariableop_19_adam_v_dense_11_biasIdentity_19:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:┤
AssignVariableOp_20AssignVariableOpassignvariableop_20_total_1Identity_20:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:┤
AssignVariableOp_21AssignVariableOpassignvariableop_21_count_1Identity_21:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:▓
AssignVariableOp_22AssignVariableOpassignvariableop_22_totalIdentity_22:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:▓
AssignVariableOp_23AssignVariableOpassignvariableop_23_countIdentity_23:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0Y
NoOpNoOp"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 ▀
Identity_24Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: W
Identity_25IdentityIdentity_24:output:0^NoOp_1*
T0*
_output_shapes
: ╠
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
_user_specified_namefile_prefix"є
L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*▒
serving_defaultЮ
;
input_40
serving_default_input_4:0         B

tf.stack_34
StatefulPartitionedCall:0         tensorflow/serving/predict:­┤
Й
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
Ъ
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
╩
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
Н
&trace_0
'trace_1
(trace_2
)trace_32Ж
(__inference_model_3_layer_call_fn_172621
(__inference_model_3_layer_call_fn_172809
(__inference_model_3_layer_call_fn_172826
(__inference_model_3_layer_call_fn_172705┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 z&trace_0z'trace_1z(trace_2z)trace_3
┴
*trace_0
+trace_1
,trace_2
-trace_32о
C__inference_model_3_layer_call_and_return_conditional_losses_172873
C__inference_model_3_layer_call_and_return_conditional_losses_172920
C__inference_model_3_layer_call_and_return_conditional_losses_172738
C__inference_model_3_layer_call_and_return_conditional_losses_172771┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 z*trace_0z+trace_1z,trace_2z-trace_3
╠B╔
!__inference__wrapped_model_172358input_4"ў
Љ▓Ї
FullArgSpec
argsџ 
varargsjargs
varkwjkwargs
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
ю
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
╗
6	variables
7trainable_variables
8regularization_losses
9	keras_api
:__call__
*;&call_and_return_all_conditional_losses

kernel
bias"
_tf_keras_layer
╗
<	variables
=trainable_variables
>regularization_losses
?	keras_api
@__call__
*A&call_and_return_all_conditional_losses

kernel
bias"
_tf_keras_layer
╗
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
Г
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
ж
Mtrace_0
Ntrace_1
Otrace_2
Ptrace_32■
-__inference_sequential_3_layer_call_fn_172431
-__inference_sequential_3_layer_call_fn_172967
-__inference_sequential_3_layer_call_fn_172984
-__inference_sequential_3_layer_call_fn_172531┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 zMtrace_0zNtrace_1zOtrace_2zPtrace_3
Н
Qtrace_0
Rtrace_1
Strace_2
Ttrace_32Ж
H__inference_sequential_3_layer_call_and_return_conditional_losses_173008
H__inference_sequential_3_layer_call_and_return_conditional_losses_173032
H__inference_sequential_3_layer_call_and_return_conditional_losses_172550
H__inference_sequential_3_layer_call_and_return_conditional_losses_172569┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 zQtrace_0zRtrace_1zStrace_2zTtrace_3
"
_generic_user_object
 :22dense_9/kernel
:22dense_9/bias
!:222dense_10/kernel
:22dense_10/bias
!:22dense_11/kernel
:2dense_11/bias
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
ЩBэ
(__inference_model_3_layer_call_fn_172621input_4"┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
щBШ
(__inference_model_3_layer_call_fn_172809inputs"┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
щBШ
(__inference_model_3_layer_call_fn_172826inputs"┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
ЩBэ
(__inference_model_3_layer_call_fn_172705input_4"┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
ћBЉ
C__inference_model_3_layer_call_and_return_conditional_losses_172873inputs"┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
ћBЉ
C__inference_model_3_layer_call_and_return_conditional_losses_172920inputs"┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
ЋBњ
C__inference_model_3_layer_call_and_return_conditional_losses_172738input_4"┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
ЋBњ
C__inference_model_3_layer_call_and_return_conditional_losses_172771input_4"┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
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
╣
ctrace_0
dtrace_1
etrace_2
ftrace_3
gtrace_4
htrace_52џ
#__inference__update_step_xla_172925
#__inference__update_step_xla_172930
#__inference__update_step_xla_172935
#__inference__update_step_xla_172940
#__inference__update_step_xla_172945
#__inference__update_step_xla_172950╣
«▓ф
FullArgSpec2
args*џ'
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

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 0zctrace_0zdtrace_1zetrace_2zftrace_3zgtrace_4zhtrace_5
╦B╚
$__inference_signature_wrapper_172792input_4"ћ
Ї▓Ѕ
FullArgSpec
argsџ 
varargs
 
varkwjkwargs
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
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
Г
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
В
ntrace_02¤
(__inference_dense_9_layer_call_fn_173041б
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 zntrace_0
Є
otrace_02Ж
C__inference_dense_9_layer_call_and_return_conditional_losses_173052б
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
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
Г
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
ь
utrace_02л
)__inference_dense_10_layer_call_fn_173061б
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 zutrace_0
ѕ
vtrace_02в
D__inference_dense_10_layer_call_and_return_conditional_losses_173072б
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
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
Г
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
ь
|trace_02л
)__inference_dense_11_layer_call_fn_173081б
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 z|trace_0
ѕ
}trace_02в
D__inference_dense_11_layer_call_and_return_conditional_losses_173091б
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
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
ЁBѓ
-__inference_sequential_3_layer_call_fn_172431dense_9_input"┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
■Bч
-__inference_sequential_3_layer_call_fn_172967inputs"┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
■Bч
-__inference_sequential_3_layer_call_fn_172984inputs"┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
ЁBѓ
-__inference_sequential_3_layer_call_fn_172531dense_9_input"┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
ЎBќ
H__inference_sequential_3_layer_call_and_return_conditional_losses_173008inputs"┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
ЎBќ
H__inference_sequential_3_layer_call_and_return_conditional_losses_173032inputs"┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
аBЮ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172550dense_9_input"┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
аBЮ
H__inference_sequential_3_layer_call_and_return_conditional_losses_172569dense_9_input"┐
Х▓▓
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
P
~	variables
	keras_api

ђtotal

Ђcount"
_tf_keras_metric
c
ѓ	variables
Ѓ	keras_api

ёtotal

Ёcount
є
_fn_kwargs"
_tf_keras_metric
%:#22Adam/m/dense_9/kernel
%:#22Adam/v/dense_9/kernel
:22Adam/m/dense_9/bias
:22Adam/v/dense_9/bias
&:$222Adam/m/dense_10/kernel
&:$222Adam/v/dense_10/kernel
 :22Adam/m/dense_10/bias
 :22Adam/v/dense_10/bias
&:$22Adam/m/dense_11/kernel
&:$22Adam/v/dense_11/kernel
 :2Adam/m/dense_11/bias
 :2Adam/v/dense_11/bias
ЭBш
#__inference__update_step_xla_172925gradientvariable"и
«▓ф
FullArgSpec2
args*џ'
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

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
ЭBш
#__inference__update_step_xla_172930gradientvariable"и
«▓ф
FullArgSpec2
args*џ'
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

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
ЭBш
#__inference__update_step_xla_172935gradientvariable"и
«▓ф
FullArgSpec2
args*џ'
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

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
ЭBш
#__inference__update_step_xla_172940gradientvariable"и
«▓ф
FullArgSpec2
args*џ'
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

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
ЭBш
#__inference__update_step_xla_172945gradientvariable"и
«▓ф
FullArgSpec2
args*џ'
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

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
ЭBш
#__inference__update_step_xla_172950gradientvariable"и
«▓ф
FullArgSpec2
args*џ'
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

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
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
▄B┘
(__inference_dense_9_layer_call_fn_173041inputs"б
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
эBЗ
C__inference_dense_9_layer_call_and_return_conditional_losses_173052inputs"б
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
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
ПB┌
)__inference_dense_10_layer_call_fn_173061inputs"б
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
ЭBш
D__inference_dense_10_layer_call_and_return_conditional_losses_173072inputs"б
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
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
ПB┌
)__inference_dense_11_layer_call_fn_173081inputs"б
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
ЭBш
D__inference_dense_11_layer_call_and_return_conditional_losses_173091inputs"б
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
0
ђ0
Ђ1"
trackable_list_wrapper
-
~	variables"
_generic_user_object
:  (2total
:  (2count
0
ё0
Ё1"
trackable_list_wrapper
.
ѓ	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapperЋ
#__inference__update_step_xla_172925nhбe
^б[
і
gradient2
4њ1	б
Щ2
ђ
p
` VariableSpec 
`ђ»ћ┬█▄?
ф "
 Ї
#__inference__update_step_xla_172930f`б]
VбS
і
gradient2
0њ-	б
Щ2
ђ
p
` VariableSpec 
`апћ┬█▄?
ф "
 Ћ
#__inference__update_step_xla_172935nhбe
^б[
і
gradient22
4њ1	б
Щ22
ђ
p
` VariableSpec 
`аСЋ┬█▄?
ф "
 Ї
#__inference__update_step_xla_172940f`б]
VбS
і
gradient2
0њ-	б
Щ2
ђ
p
` VariableSpec 
`ђсЋ┬█▄?
ф "
 Ћ
#__inference__update_step_xla_172945nhбe
^б[
і
gradient2
4њ1	б
Щ2
ђ
p
` VariableSpec 
`ђ├ћ┬█▄?
ф "
 Ї
#__inference__update_step_xla_172950f`б]
VбS
і
gradient
0њ-	б
Щ
ђ
p
` VariableSpec 
`аПћ┬█▄?
ф "
 ю
!__inference__wrapped_model_172358w 0б-
&б#
!і
input_4         
ф ";ф8
6

tf.stack_3(і%

tf_stack_3         Ф
D__inference_dense_10_layer_call_and_return_conditional_losses_173072c/б,
%б"
 і
inputs         2
ф ",б)
"і
tensor_0         2
џ Ё
)__inference_dense_10_layer_call_fn_173061X/б,
%б"
 і
inputs         2
ф "!і
unknown         2Ф
D__inference_dense_11_layer_call_and_return_conditional_losses_173091c /б,
%б"
 і
inputs         2
ф ",б)
"і
tensor_0         
џ Ё
)__inference_dense_11_layer_call_fn_173081X /б,
%б"
 і
inputs         2
ф "!і
unknown         ф
C__inference_dense_9_layer_call_and_return_conditional_losses_173052c/б,
%б"
 і
inputs         
ф ",б)
"і
tensor_0         2
џ ё
(__inference_dense_9_layer_call_fn_173041X/б,
%б"
 і
inputs         
ф "!і
unknown         2╗
C__inference_model_3_layer_call_and_return_conditional_losses_172738t 8б5
.б+
!і
input_4         
p 

 
ф "0б-
&і#
tensor_0         
џ ╗
C__inference_model_3_layer_call_and_return_conditional_losses_172771t 8б5
.б+
!і
input_4         
p

 
ф "0б-
&і#
tensor_0         
џ ║
C__inference_model_3_layer_call_and_return_conditional_losses_172873s 7б4
-б*
 і
inputs         
p 

 
ф "0б-
&і#
tensor_0         
џ ║
C__inference_model_3_layer_call_and_return_conditional_losses_172920s 7б4
-б*
 і
inputs         
p

 
ф "0б-
&і#
tensor_0         
џ Ћ
(__inference_model_3_layer_call_fn_172621i 8б5
.б+
!і
input_4         
p 

 
ф "%і"
unknown         Ћ
(__inference_model_3_layer_call_fn_172705i 8б5
.б+
!і
input_4         
p

 
ф "%і"
unknown         ћ
(__inference_model_3_layer_call_fn_172809h 7б4
-б*
 і
inputs         
p 

 
ф "%і"
unknown         ћ
(__inference_model_3_layer_call_fn_172826h 7б4
-б*
 і
inputs         
p

 
ф "%і"
unknown         ┬
H__inference_sequential_3_layer_call_and_return_conditional_losses_172550v >б;
4б1
'і$
dense_9_input         
p 

 
ф ",б)
"і
tensor_0         
џ ┬
H__inference_sequential_3_layer_call_and_return_conditional_losses_172569v >б;
4б1
'і$
dense_9_input         
p

 
ф ",б)
"і
tensor_0         
џ ╗
H__inference_sequential_3_layer_call_and_return_conditional_losses_173008o 7б4
-б*
 і
inputs         
p 

 
ф ",б)
"і
tensor_0         
џ ╗
H__inference_sequential_3_layer_call_and_return_conditional_losses_173032o 7б4
-б*
 і
inputs         
p

 
ф ",б)
"і
tensor_0         
џ ю
-__inference_sequential_3_layer_call_fn_172431k >б;
4б1
'і$
dense_9_input         
p 

 
ф "!і
unknown         ю
-__inference_sequential_3_layer_call_fn_172531k >б;
4б1
'і$
dense_9_input         
p

 
ф "!і
unknown         Ћ
-__inference_sequential_3_layer_call_fn_172967d 7б4
-б*
 і
inputs         
p 

 
ф "!і
unknown         Ћ
-__inference_sequential_3_layer_call_fn_172984d 7б4
-б*
 і
inputs         
p

 
ф "!і
unknown         Ф
$__inference_signature_wrapper_172792ѓ ;б8
б 
1ф.
,
input_4!і
input_4         ";ф8
6

tf.stack_3(і%

tf_stack_3         