Иє
НЁ
^
AssignVariableOp
resource
value"dtype"
dtypetype"
validate_shapebool( И
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
Ж
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool("
allow_missing_filesbool( И
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
dtypetypeИ
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
list(type)(0И
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0И
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
0
Sigmoid
x"T
y"T"
Ttype:

2
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
executor_typestring Ии
@
StaticRegexFullMatch	
input

output
"
patternstring
ў
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
Ц
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 И"serve*2.11.02unknown8▒к
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
А
Adam/v/dense_14/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/v/dense_14/bias
y
(Adam/v/dense_14/bias/Read/ReadVariableOpReadVariableOpAdam/v/dense_14/bias*
_output_shapes
:*
dtype0
А
Adam/m/dense_14/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/m/dense_14/bias
y
(Adam/m/dense_14/bias/Read/ReadVariableOpReadVariableOpAdam/m/dense_14/bias*
_output_shapes
:*
dtype0
И
Adam/v/dense_14/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*'
shared_nameAdam/v/dense_14/kernel
Б
*Adam/v/dense_14/kernel/Read/ReadVariableOpReadVariableOpAdam/v/dense_14/kernel*
_output_shapes

:2*
dtype0
И
Adam/m/dense_14/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*'
shared_nameAdam/m/dense_14/kernel
Б
*Adam/m/dense_14/kernel/Read/ReadVariableOpReadVariableOpAdam/m/dense_14/kernel*
_output_shapes

:2*
dtype0
А
Adam/v/dense_13/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*%
shared_nameAdam/v/dense_13/bias
y
(Adam/v/dense_13/bias/Read/ReadVariableOpReadVariableOpAdam/v/dense_13/bias*
_output_shapes
:2*
dtype0
А
Adam/m/dense_13/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*%
shared_nameAdam/m/dense_13/bias
y
(Adam/m/dense_13/bias/Read/ReadVariableOpReadVariableOpAdam/m/dense_13/bias*
_output_shapes
:2*
dtype0
И
Adam/v/dense_13/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*'
shared_nameAdam/v/dense_13/kernel
Б
*Adam/v/dense_13/kernel/Read/ReadVariableOpReadVariableOpAdam/v/dense_13/kernel*
_output_shapes

:22*
dtype0
И
Adam/m/dense_13/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*'
shared_nameAdam/m/dense_13/kernel
Б
*Adam/m/dense_13/kernel/Read/ReadVariableOpReadVariableOpAdam/m/dense_13/kernel*
_output_shapes

:22*
dtype0
А
Adam/v/dense_12/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*%
shared_nameAdam/v/dense_12/bias
y
(Adam/v/dense_12/bias/Read/ReadVariableOpReadVariableOpAdam/v/dense_12/bias*
_output_shapes
:2*
dtype0
А
Adam/m/dense_12/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*%
shared_nameAdam/m/dense_12/bias
y
(Adam/m/dense_12/bias/Read/ReadVariableOpReadVariableOpAdam/m/dense_12/bias*
_output_shapes
:2*
dtype0
И
Adam/v/dense_12/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:	2*'
shared_nameAdam/v/dense_12/kernel
Б
*Adam/v/dense_12/kernel/Read/ReadVariableOpReadVariableOpAdam/v/dense_12/kernel*
_output_shapes

:	2*
dtype0
И
Adam/m/dense_12/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:	2*'
shared_nameAdam/m/dense_12/kernel
Б
*Adam/m/dense_12/kernel/Read/ReadVariableOpReadVariableOpAdam/m/dense_12/kernel*
_output_shapes

:	2*
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
dense_14/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_namedense_14/bias
k
!dense_14/bias/Read/ReadVariableOpReadVariableOpdense_14/bias*
_output_shapes
:*
dtype0
z
dense_14/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2* 
shared_namedense_14/kernel
s
#dense_14/kernel/Read/ReadVariableOpReadVariableOpdense_14/kernel*
_output_shapes

:2*
dtype0
r
dense_13/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*
shared_namedense_13/bias
k
!dense_13/bias/Read/ReadVariableOpReadVariableOpdense_13/bias*
_output_shapes
:2*
dtype0
z
dense_13/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22* 
shared_namedense_13/kernel
s
#dense_13/kernel/Read/ReadVariableOpReadVariableOpdense_13/kernel*
_output_shapes

:22*
dtype0
r
dense_12/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*
shared_namedense_12/bias
k
!dense_12/bias/Read/ReadVariableOpReadVariableOpdense_12/bias*
_output_shapes
:2*
dtype0
z
dense_12/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:	2* 
shared_namedense_12/kernel
s
#dense_12/kernel/Read/ReadVariableOpReadVariableOpdense_12/kernel*
_output_shapes

:	2*
dtype0
z
serving_default_input_5Placeholder*'
_output_shapes
:         *
dtype0*
shape:         
и
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_5dense_12/kerneldense_12/biasdense_13/kerneldense_13/biasdense_14/kerneldense_14/bias*
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
GPU2 *0J 8В *.
f)R'
%__inference_signature_wrapper_7019352

NoOpNoOp
╡2
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*Ё1
valueц1Bу1 B▄1
з
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
Е
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
Б
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
ж
6	variables
7trainable_variables
8regularization_losses
9	keras_api
:__call__
*;&call_and_return_all_conditional_losses

kernel
bias*
ж
<	variables
=trainable_variables
>regularization_losses
?	keras_api
@__call__
*A&call_and_return_all_conditional_losses

kernel
bias*
ж
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
У
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
OI
VARIABLE_VALUEdense_12/kernel&variables/0/.ATTRIBUTES/VARIABLE_VALUE*
MG
VARIABLE_VALUEdense_12/bias&variables/1/.ATTRIBUTES/VARIABLE_VALUE*
OI
VARIABLE_VALUEdense_13/kernel&variables/2/.ATTRIBUTES/VARIABLE_VALUE*
MG
VARIABLE_VALUEdense_13/bias&variables/3/.ATTRIBUTES/VARIABLE_VALUE*
OI
VARIABLE_VALUEdense_14/kernel&variables/4/.ATTRIBUTES/VARIABLE_VALUE*
MG
VARIABLE_VALUEdense_14/bias&variables/5/.ATTRIBUTES/VARIABLE_VALUE*
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
У
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
У
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
У
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

Аtotal

Бcount*
M
В	variables
Г	keras_api

Дtotal

Еcount
Ж
_fn_kwargs*
a[
VARIABLE_VALUEAdam/m/dense_12/kernel1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/v/dense_12/kernel1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/m/dense_12/bias1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/v/dense_12/bias1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/m/dense_13/kernel1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/v/dense_13/kernel1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/m/dense_13/bias1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/v/dense_13/bias1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/m/dense_14/kernel1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUE*
b\
VARIABLE_VALUEAdam/v/dense_14/kernel2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUEAdam/m/dense_14/bias2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUEAdam/v/dense_14/bias2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUE*
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
А0
Б1*

~	variables*
UO
VARIABLE_VALUEtotal_14keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEcount_14keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE*

Д0
Е1*

В	variables*
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
╟	
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename#dense_12/kernel/Read/ReadVariableOp!dense_12/bias/Read/ReadVariableOp#dense_13/kernel/Read/ReadVariableOp!dense_13/bias/Read/ReadVariableOp#dense_14/kernel/Read/ReadVariableOp!dense_14/bias/Read/ReadVariableOpiteration/Read/ReadVariableOp!learning_rate/Read/ReadVariableOp*Adam/m/dense_12/kernel/Read/ReadVariableOp*Adam/v/dense_12/kernel/Read/ReadVariableOp(Adam/m/dense_12/bias/Read/ReadVariableOp(Adam/v/dense_12/bias/Read/ReadVariableOp*Adam/m/dense_13/kernel/Read/ReadVariableOp*Adam/v/dense_13/kernel/Read/ReadVariableOp(Adam/m/dense_13/bias/Read/ReadVariableOp(Adam/v/dense_13/bias/Read/ReadVariableOp*Adam/m/dense_14/kernel/Read/ReadVariableOp*Adam/v/dense_14/kernel/Read/ReadVariableOp(Adam/m/dense_14/bias/Read/ReadVariableOp(Adam/v/dense_14/bias/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOpConst*%
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
GPU2 *0J 8В *)
f$R"
 __inference__traced_save_7019753
т
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenamedense_12/kerneldense_12/biasdense_13/kerneldense_13/biasdense_14/kerneldense_14/bias	iterationlearning_rateAdam/m/dense_12/kernelAdam/v/dense_12/kernelAdam/m/dense_12/biasAdam/v/dense_12/biasAdam/m/dense_13/kernelAdam/v/dense_13/kernelAdam/m/dense_13/biasAdam/v/dense_13/biasAdam/m/dense_14/kernelAdam/v/dense_14/kernelAdam/m/dense_14/biasAdam/v/dense_14/biastotal_1count_1totalcount*$
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
GPU2 *0J 8В *,
f'R%
#__inference__traced_restore_7019835▀╢
н
L
$__inference__update_step_xla_7019504
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
е
К
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019598

inputs9
'dense_12_matmul_readvariableop_resource:	26
(dense_12_biasadd_readvariableop_resource:29
'dense_13_matmul_readvariableop_resource:226
(dense_13_biasadd_readvariableop_resource:29
'dense_14_matmul_readvariableop_resource:26
(dense_14_biasadd_readvariableop_resource:
identityИвdense_12/BiasAdd/ReadVariableOpвdense_12/MatMul/ReadVariableOpвdense_13/BiasAdd/ReadVariableOpвdense_13/MatMul/ReadVariableOpвdense_14/BiasAdd/ReadVariableOpвdense_14/MatMul/ReadVariableOpЖ
dense_12/MatMul/ReadVariableOpReadVariableOp'dense_12_matmul_readvariableop_resource*
_output_shapes

:	2*
dtype0{
dense_12/MatMulMatMulinputs&dense_12/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2Д
dense_12/BiasAdd/ReadVariableOpReadVariableOp(dense_12_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0С
dense_12/BiasAddBiasAdddense_12/MatMul:product:0'dense_12/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2b
dense_12/ReluReludense_12/BiasAdd:output:0*
T0*'
_output_shapes
:         2Ж
dense_13/MatMul/ReadVariableOpReadVariableOp'dense_13_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Р
dense_13/MatMulMatMuldense_12/Relu:activations:0&dense_13/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2Д
dense_13/BiasAdd/ReadVariableOpReadVariableOp(dense_13_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0С
dense_13/BiasAddBiasAdddense_13/MatMul:product:0'dense_13/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2b
dense_13/ReluReludense_13/BiasAdd:output:0*
T0*'
_output_shapes
:         2Ж
dense_14/MatMul/ReadVariableOpReadVariableOp'dense_14_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0Р
dense_14/MatMulMatMuldense_13/Relu:activations:0&dense_14/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         Д
dense_14/BiasAdd/ReadVariableOpReadVariableOp(dense_14_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0С
dense_14/BiasAddBiasAdddense_14/MatMul:product:0'dense_14/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         h
dense_14/SigmoidSigmoiddense_14/BiasAdd:output:0*
T0*'
_output_shapes
:         c
IdentityIdentitydense_14/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:         П
NoOpNoOp ^dense_12/BiasAdd/ReadVariableOp^dense_12/MatMul/ReadVariableOp ^dense_13/BiasAdd/ReadVariableOp^dense_13/MatMul/ReadVariableOp ^dense_14/BiasAdd/ReadVariableOp^dense_14/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         	: : : : : : 2B
dense_12/BiasAdd/ReadVariableOpdense_12/BiasAdd/ReadVariableOp2@
dense_12/MatMul/ReadVariableOpdense_12/MatMul/ReadVariableOp2B
dense_13/BiasAdd/ReadVariableOpdense_13/BiasAdd/ReadVariableOp2@
dense_13/MatMul/ReadVariableOpdense_13/MatMul/ReadVariableOp2B
dense_14/BiasAdd/ReadVariableOpdense_14/BiasAdd/ReadVariableOp2@
dense_14/MatMul/ReadVariableOpdense_14/MatMul/ReadVariableOp:O K
'
_output_shapes
:         	
 
_user_specified_nameinputs
·
Г
)__inference_model_4_layer_call_fn_7019181
input_5
unknown:	2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityИвStatefulPartitionedCallЧ
StatefulPartitionedCallStatefulPartitionedCallinput_5unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
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
GPU2 *0J 8В *M
fHRF
D__inference_model_4_layer_call_and_return_conditional_losses_7019166s
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
:         : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:         
!
_user_specified_name	input_5
╣
P
$__inference__update_step_xla_7019489
gradient
variable:	2*
_XlaMustCompile(*(
_construction_contextkEagerRuntime*
_input_shapes
:	2: *
	_noinline(:H D

_output_shapes

:	2
"
_user_specified_name
gradient:($
"
_user_specified_name
variable
·
Г
)__inference_model_4_layer_call_fn_7019265
input_5
unknown:	2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityИвStatefulPartitionedCallЧ
StatefulPartitionedCallStatefulPartitionedCallinput_5unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
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
GPU2 *0J 8В *M
fHRF
D__inference_model_4_layer_call_and_return_conditional_losses_7019233s
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
:         : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:         
!
_user_specified_name	input_5
┐
Ш
D__inference_model_4_layer_call_and_return_conditional_losses_7019331
input_5&
sequential_4_7019309:	2"
sequential_4_7019311:2&
sequential_4_7019313:22"
sequential_4_7019315:2&
sequential_4_7019317:2"
sequential_4_7019319:
identityИв$sequential_4/StatefulPartitionedCallв&sequential_4/StatefulPartitionedCall_1
.tf.__operators__.getitem_8/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"    	   Б
0tf.__operators__.getitem_8/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        Б
0tf.__operators__.getitem_8/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ╤
(tf.__operators__.getitem_8/strided_sliceStridedSliceinput_57tf.__operators__.getitem_8/strided_slice/stack:output:09tf.__operators__.getitem_8/strided_slice/stack_1:output:09tf.__operators__.getitem_8/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         	*

begin_mask*
end_mask
.tf.__operators__.getitem_7/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        Б
0tf.__operators__.getitem_7/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"    	   Б
0tf.__operators__.getitem_7/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ╤
(tf.__operators__.getitem_7/strided_sliceStridedSliceinput_57tf.__operators__.getitem_7/strided_slice/stack:output:09tf.__operators__.getitem_7/strided_slice/stack_1:output:09tf.__operators__.getitem_7/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         	*

begin_mask*
end_maskУ
$sequential_4/StatefulPartitionedCallStatefulPartitionedCall1tf.__operators__.getitem_7/strided_slice:output:0sequential_4_7019309sequential_4_7019311sequential_4_7019313sequential_4_7019315sequential_4_7019317sequential_4_7019319*
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
GPU2 *0J 8В *R
fMRK
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019059Х
&sequential_4/StatefulPartitionedCall_1StatefulPartitionedCall1tf.__operators__.getitem_8/strided_slice:output:0sequential_4_7019309sequential_4_7019311sequential_4_7019313sequential_4_7019315sequential_4_7019317sequential_4_7019319*
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
GPU2 *0J 8В *R
fMRK
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019059├
tf.stack_3/stackPack-sequential_4/StatefulPartitionedCall:output:0/sequential_4/StatefulPartitionedCall_1:output:0*
N*
T0*+
_output_shapes
:         *

axisl
IdentityIdentitytf.stack_3/stack:output:0^NoOp*
T0*+
_output_shapes
:         Ц
NoOpNoOp%^sequential_4/StatefulPartitionedCall'^sequential_4/StatefulPartitionedCall_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2L
$sequential_4/StatefulPartitionedCall$sequential_4/StatefulPartitionedCall2P
&sequential_4/StatefulPartitionedCall_1&sequential_4/StatefulPartitionedCall_1:P L
'
_output_shapes
:         
!
_user_specified_name	input_5
╘
 
%__inference_signature_wrapper_7019352
input_5
unknown:	2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityИвStatefulPartitionedCallї
StatefulPartitionedCallStatefulPartitionedCallinput_5unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
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
GPU2 *0J 8В *+
f&R$
"__inference__wrapped_model_7018917s
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
:         : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:         
!
_user_specified_name	input_5
С	
П
.__inference_sequential_4_layer_call_fn_7018991
dense_12_input
unknown:	2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityИвStatefulPartitionedCallЯ
StatefulPartitionedCallStatefulPartitionedCalldense_12_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
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
GPU2 *0J 8В *R
fMRK
I__inference_sequential_4_layer_call_and_return_conditional_losses_7018976o
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
:         	: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
'
_output_shapes
:         	
(
_user_specified_namedense_12_input
Ь

Ў
E__inference_dense_12_layer_call_and_return_conditional_losses_7019618

inputs0
matmul_readvariableop_resource:	2-
biasadd_readvariableop_resource:2
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:	2*
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
:         	: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:         	
 
_user_specified_nameinputs
ЖM
─
D__inference_model_4_layer_call_and_return_conditional_losses_7019484

inputsF
4sequential_4_dense_12_matmul_readvariableop_resource:	2C
5sequential_4_dense_12_biasadd_readvariableop_resource:2F
4sequential_4_dense_13_matmul_readvariableop_resource:22C
5sequential_4_dense_13_biasadd_readvariableop_resource:2F
4sequential_4_dense_14_matmul_readvariableop_resource:2C
5sequential_4_dense_14_biasadd_readvariableop_resource:
identityИв,sequential_4/dense_12/BiasAdd/ReadVariableOpв.sequential_4/dense_12/BiasAdd_1/ReadVariableOpв+sequential_4/dense_12/MatMul/ReadVariableOpв-sequential_4/dense_12/MatMul_1/ReadVariableOpв,sequential_4/dense_13/BiasAdd/ReadVariableOpв.sequential_4/dense_13/BiasAdd_1/ReadVariableOpв+sequential_4/dense_13/MatMul/ReadVariableOpв-sequential_4/dense_13/MatMul_1/ReadVariableOpв,sequential_4/dense_14/BiasAdd/ReadVariableOpв.sequential_4/dense_14/BiasAdd_1/ReadVariableOpв+sequential_4/dense_14/MatMul/ReadVariableOpв-sequential_4/dense_14/MatMul_1/ReadVariableOp
.tf.__operators__.getitem_8/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"    	   Б
0tf.__operators__.getitem_8/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        Б
0tf.__operators__.getitem_8/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ╨
(tf.__operators__.getitem_8/strided_sliceStridedSliceinputs7tf.__operators__.getitem_8/strided_slice/stack:output:09tf.__operators__.getitem_8/strided_slice/stack_1:output:09tf.__operators__.getitem_8/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         	*

begin_mask*
end_mask
.tf.__operators__.getitem_7/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        Б
0tf.__operators__.getitem_7/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"    	   Б
0tf.__operators__.getitem_7/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ╨
(tf.__operators__.getitem_7/strided_sliceStridedSliceinputs7tf.__operators__.getitem_7/strided_slice/stack:output:09tf.__operators__.getitem_7/strided_slice/stack_1:output:09tf.__operators__.getitem_7/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         	*

begin_mask*
end_maskа
+sequential_4/dense_12/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_12_matmul_readvariableop_resource*
_output_shapes

:	2*
dtype0└
sequential_4/dense_12/MatMulMatMul1tf.__operators__.getitem_7/strided_slice:output:03sequential_4/dense_12/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2Ю
,sequential_4/dense_12/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_12_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0╕
sequential_4/dense_12/BiasAddBiasAdd&sequential_4/dense_12/MatMul:product:04sequential_4/dense_12/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2|
sequential_4/dense_12/ReluRelu&sequential_4/dense_12/BiasAdd:output:0*
T0*'
_output_shapes
:         2а
+sequential_4/dense_13/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_13_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0╖
sequential_4/dense_13/MatMulMatMul(sequential_4/dense_12/Relu:activations:03sequential_4/dense_13/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2Ю
,sequential_4/dense_13/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_13_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0╕
sequential_4/dense_13/BiasAddBiasAdd&sequential_4/dense_13/MatMul:product:04sequential_4/dense_13/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2|
sequential_4/dense_13/ReluRelu&sequential_4/dense_13/BiasAdd:output:0*
T0*'
_output_shapes
:         2а
+sequential_4/dense_14/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_14_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0╖
sequential_4/dense_14/MatMulMatMul(sequential_4/dense_13/Relu:activations:03sequential_4/dense_14/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         Ю
,sequential_4/dense_14/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_14_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0╕
sequential_4/dense_14/BiasAddBiasAdd&sequential_4/dense_14/MatMul:product:04sequential_4/dense_14/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         В
sequential_4/dense_14/SigmoidSigmoid&sequential_4/dense_14/BiasAdd:output:0*
T0*'
_output_shapes
:         в
-sequential_4/dense_12/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_12_matmul_readvariableop_resource*
_output_shapes

:	2*
dtype0─
sequential_4/dense_12/MatMul_1MatMul1tf.__operators__.getitem_8/strided_slice:output:05sequential_4/dense_12/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2а
.sequential_4/dense_12/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_12_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0╛
sequential_4/dense_12/BiasAdd_1BiasAdd(sequential_4/dense_12/MatMul_1:product:06sequential_4/dense_12/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2А
sequential_4/dense_12/Relu_1Relu(sequential_4/dense_12/BiasAdd_1:output:0*
T0*'
_output_shapes
:         2в
-sequential_4/dense_13/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_13_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0╜
sequential_4/dense_13/MatMul_1MatMul*sequential_4/dense_12/Relu_1:activations:05sequential_4/dense_13/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2а
.sequential_4/dense_13/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_13_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0╛
sequential_4/dense_13/BiasAdd_1BiasAdd(sequential_4/dense_13/MatMul_1:product:06sequential_4/dense_13/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2А
sequential_4/dense_13/Relu_1Relu(sequential_4/dense_13/BiasAdd_1:output:0*
T0*'
_output_shapes
:         2в
-sequential_4/dense_14/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_14_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0╜
sequential_4/dense_14/MatMul_1MatMul*sequential_4/dense_13/Relu_1:activations:05sequential_4/dense_14/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         а
.sequential_4/dense_14/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_14_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0╛
sequential_4/dense_14/BiasAdd_1BiasAdd(sequential_4/dense_14/MatMul_1:product:06sequential_4/dense_14/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         Ж
sequential_4/dense_14/Sigmoid_1Sigmoid(sequential_4/dense_14/BiasAdd_1:output:0*
T0*'
_output_shapes
:         л
tf.stack_3/stackPack!sequential_4/dense_14/Sigmoid:y:0#sequential_4/dense_14/Sigmoid_1:y:0*
N*
T0*+
_output_shapes
:         *

axisl
IdentityIdentitytf.stack_3/stack:output:0^NoOp*
T0*+
_output_shapes
:         А
NoOpNoOp-^sequential_4/dense_12/BiasAdd/ReadVariableOp/^sequential_4/dense_12/BiasAdd_1/ReadVariableOp,^sequential_4/dense_12/MatMul/ReadVariableOp.^sequential_4/dense_12/MatMul_1/ReadVariableOp-^sequential_4/dense_13/BiasAdd/ReadVariableOp/^sequential_4/dense_13/BiasAdd_1/ReadVariableOp,^sequential_4/dense_13/MatMul/ReadVariableOp.^sequential_4/dense_13/MatMul_1/ReadVariableOp-^sequential_4/dense_14/BiasAdd/ReadVariableOp/^sequential_4/dense_14/BiasAdd_1/ReadVariableOp,^sequential_4/dense_14/MatMul/ReadVariableOp.^sequential_4/dense_14/MatMul_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2\
,sequential_4/dense_12/BiasAdd/ReadVariableOp,sequential_4/dense_12/BiasAdd/ReadVariableOp2`
.sequential_4/dense_12/BiasAdd_1/ReadVariableOp.sequential_4/dense_12/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_12/MatMul/ReadVariableOp+sequential_4/dense_12/MatMul/ReadVariableOp2^
-sequential_4/dense_12/MatMul_1/ReadVariableOp-sequential_4/dense_12/MatMul_1/ReadVariableOp2\
,sequential_4/dense_13/BiasAdd/ReadVariableOp,sequential_4/dense_13/BiasAdd/ReadVariableOp2`
.sequential_4/dense_13/BiasAdd_1/ReadVariableOp.sequential_4/dense_13/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_13/MatMul/ReadVariableOp+sequential_4/dense_13/MatMul/ReadVariableOp2^
-sequential_4/dense_13/MatMul_1/ReadVariableOp-sequential_4/dense_13/MatMul_1/ReadVariableOp2\
,sequential_4/dense_14/BiasAdd/ReadVariableOp,sequential_4/dense_14/BiasAdd/ReadVariableOp2`
.sequential_4/dense_14/BiasAdd_1/ReadVariableOp.sequential_4/dense_14/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_14/MatMul/ReadVariableOp+sequential_4/dense_14/MatMul/ReadVariableOp2^
-sequential_4/dense_14/MatMul_1/ReadVariableOp-sequential_4/dense_14/MatMul_1/ReadVariableOp:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
ў
В
)__inference_model_4_layer_call_fn_7019386

inputs
unknown:	2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityИвStatefulPartitionedCallЦ
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
GPU2 *0J 8В *M
fHRF
D__inference_model_4_layer_call_and_return_conditional_losses_7019233s
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
:         : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
Ь

Ў
E__inference_dense_13_layer_call_and_return_conditional_losses_7018952

inputs0
matmul_readvariableop_resource:22-
biasadd_readvariableop_resource:2
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpt
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
Ь

Ў
E__inference_dense_12_layer_call_and_return_conditional_losses_7018935

inputs0
matmul_readvariableop_resource:	2-
biasadd_readvariableop_resource:2
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:	2*
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
:         	: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:         	
 
_user_specified_nameinputs
ьU
│	
"__inference__wrapped_model_7018917
input_5N
<model_4_sequential_4_dense_12_matmul_readvariableop_resource:	2K
=model_4_sequential_4_dense_12_biasadd_readvariableop_resource:2N
<model_4_sequential_4_dense_13_matmul_readvariableop_resource:22K
=model_4_sequential_4_dense_13_biasadd_readvariableop_resource:2N
<model_4_sequential_4_dense_14_matmul_readvariableop_resource:2K
=model_4_sequential_4_dense_14_biasadd_readvariableop_resource:
identityИв4model_4/sequential_4/dense_12/BiasAdd/ReadVariableOpв6model_4/sequential_4/dense_12/BiasAdd_1/ReadVariableOpв3model_4/sequential_4/dense_12/MatMul/ReadVariableOpв5model_4/sequential_4/dense_12/MatMul_1/ReadVariableOpв4model_4/sequential_4/dense_13/BiasAdd/ReadVariableOpв6model_4/sequential_4/dense_13/BiasAdd_1/ReadVariableOpв3model_4/sequential_4/dense_13/MatMul/ReadVariableOpв5model_4/sequential_4/dense_13/MatMul_1/ReadVariableOpв4model_4/sequential_4/dense_14/BiasAdd/ReadVariableOpв6model_4/sequential_4/dense_14/BiasAdd_1/ReadVariableOpв3model_4/sequential_4/dense_14/MatMul/ReadVariableOpв5model_4/sequential_4/dense_14/MatMul_1/ReadVariableOpЗ
6model_4/tf.__operators__.getitem_8/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"    	   Й
8model_4/tf.__operators__.getitem_8/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        Й
8model_4/tf.__operators__.getitem_8/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ё
0model_4/tf.__operators__.getitem_8/strided_sliceStridedSliceinput_5?model_4/tf.__operators__.getitem_8/strided_slice/stack:output:0Amodel_4/tf.__operators__.getitem_8/strided_slice/stack_1:output:0Amodel_4/tf.__operators__.getitem_8/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         	*

begin_mask*
end_maskЗ
6model_4/tf.__operators__.getitem_7/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        Й
8model_4/tf.__operators__.getitem_7/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"    	   Й
8model_4/tf.__operators__.getitem_7/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ё
0model_4/tf.__operators__.getitem_7/strided_sliceStridedSliceinput_5?model_4/tf.__operators__.getitem_7/strided_slice/stack:output:0Amodel_4/tf.__operators__.getitem_7/strided_slice/stack_1:output:0Amodel_4/tf.__operators__.getitem_7/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         	*

begin_mask*
end_mask░
3model_4/sequential_4/dense_12/MatMul/ReadVariableOpReadVariableOp<model_4_sequential_4_dense_12_matmul_readvariableop_resource*
_output_shapes

:	2*
dtype0╪
$model_4/sequential_4/dense_12/MatMulMatMul9model_4/tf.__operators__.getitem_7/strided_slice:output:0;model_4/sequential_4/dense_12/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2о
4model_4/sequential_4/dense_12/BiasAdd/ReadVariableOpReadVariableOp=model_4_sequential_4_dense_12_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0╨
%model_4/sequential_4/dense_12/BiasAddBiasAdd.model_4/sequential_4/dense_12/MatMul:product:0<model_4/sequential_4/dense_12/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2М
"model_4/sequential_4/dense_12/ReluRelu.model_4/sequential_4/dense_12/BiasAdd:output:0*
T0*'
_output_shapes
:         2░
3model_4/sequential_4/dense_13/MatMul/ReadVariableOpReadVariableOp<model_4_sequential_4_dense_13_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0╧
$model_4/sequential_4/dense_13/MatMulMatMul0model_4/sequential_4/dense_12/Relu:activations:0;model_4/sequential_4/dense_13/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2о
4model_4/sequential_4/dense_13/BiasAdd/ReadVariableOpReadVariableOp=model_4_sequential_4_dense_13_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0╨
%model_4/sequential_4/dense_13/BiasAddBiasAdd.model_4/sequential_4/dense_13/MatMul:product:0<model_4/sequential_4/dense_13/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2М
"model_4/sequential_4/dense_13/ReluRelu.model_4/sequential_4/dense_13/BiasAdd:output:0*
T0*'
_output_shapes
:         2░
3model_4/sequential_4/dense_14/MatMul/ReadVariableOpReadVariableOp<model_4_sequential_4_dense_14_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0╧
$model_4/sequential_4/dense_14/MatMulMatMul0model_4/sequential_4/dense_13/Relu:activations:0;model_4/sequential_4/dense_14/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         о
4model_4/sequential_4/dense_14/BiasAdd/ReadVariableOpReadVariableOp=model_4_sequential_4_dense_14_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0╨
%model_4/sequential_4/dense_14/BiasAddBiasAdd.model_4/sequential_4/dense_14/MatMul:product:0<model_4/sequential_4/dense_14/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         Т
%model_4/sequential_4/dense_14/SigmoidSigmoid.model_4/sequential_4/dense_14/BiasAdd:output:0*
T0*'
_output_shapes
:         ▓
5model_4/sequential_4/dense_12/MatMul_1/ReadVariableOpReadVariableOp<model_4_sequential_4_dense_12_matmul_readvariableop_resource*
_output_shapes

:	2*
dtype0▄
&model_4/sequential_4/dense_12/MatMul_1MatMul9model_4/tf.__operators__.getitem_8/strided_slice:output:0=model_4/sequential_4/dense_12/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2░
6model_4/sequential_4/dense_12/BiasAdd_1/ReadVariableOpReadVariableOp=model_4_sequential_4_dense_12_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0╓
'model_4/sequential_4/dense_12/BiasAdd_1BiasAdd0model_4/sequential_4/dense_12/MatMul_1:product:0>model_4/sequential_4/dense_12/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2Р
$model_4/sequential_4/dense_12/Relu_1Relu0model_4/sequential_4/dense_12/BiasAdd_1:output:0*
T0*'
_output_shapes
:         2▓
5model_4/sequential_4/dense_13/MatMul_1/ReadVariableOpReadVariableOp<model_4_sequential_4_dense_13_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0╒
&model_4/sequential_4/dense_13/MatMul_1MatMul2model_4/sequential_4/dense_12/Relu_1:activations:0=model_4/sequential_4/dense_13/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2░
6model_4/sequential_4/dense_13/BiasAdd_1/ReadVariableOpReadVariableOp=model_4_sequential_4_dense_13_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0╓
'model_4/sequential_4/dense_13/BiasAdd_1BiasAdd0model_4/sequential_4/dense_13/MatMul_1:product:0>model_4/sequential_4/dense_13/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2Р
$model_4/sequential_4/dense_13/Relu_1Relu0model_4/sequential_4/dense_13/BiasAdd_1:output:0*
T0*'
_output_shapes
:         2▓
5model_4/sequential_4/dense_14/MatMul_1/ReadVariableOpReadVariableOp<model_4_sequential_4_dense_14_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0╒
&model_4/sequential_4/dense_14/MatMul_1MatMul2model_4/sequential_4/dense_13/Relu_1:activations:0=model_4/sequential_4/dense_14/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         ░
6model_4/sequential_4/dense_14/BiasAdd_1/ReadVariableOpReadVariableOp=model_4_sequential_4_dense_14_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0╓
'model_4/sequential_4/dense_14/BiasAdd_1BiasAdd0model_4/sequential_4/dense_14/MatMul_1:product:0>model_4/sequential_4/dense_14/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         Ц
'model_4/sequential_4/dense_14/Sigmoid_1Sigmoid0model_4/sequential_4/dense_14/BiasAdd_1:output:0*
T0*'
_output_shapes
:         ├
model_4/tf.stack_3/stackPack)model_4/sequential_4/dense_14/Sigmoid:y:0+model_4/sequential_4/dense_14/Sigmoid_1:y:0*
N*
T0*+
_output_shapes
:         *

axist
IdentityIdentity!model_4/tf.stack_3/stack:output:0^NoOp*
T0*+
_output_shapes
:         р
NoOpNoOp5^model_4/sequential_4/dense_12/BiasAdd/ReadVariableOp7^model_4/sequential_4/dense_12/BiasAdd_1/ReadVariableOp4^model_4/sequential_4/dense_12/MatMul/ReadVariableOp6^model_4/sequential_4/dense_12/MatMul_1/ReadVariableOp5^model_4/sequential_4/dense_13/BiasAdd/ReadVariableOp7^model_4/sequential_4/dense_13/BiasAdd_1/ReadVariableOp4^model_4/sequential_4/dense_13/MatMul/ReadVariableOp6^model_4/sequential_4/dense_13/MatMul_1/ReadVariableOp5^model_4/sequential_4/dense_14/BiasAdd/ReadVariableOp7^model_4/sequential_4/dense_14/BiasAdd_1/ReadVariableOp4^model_4/sequential_4/dense_14/MatMul/ReadVariableOp6^model_4/sequential_4/dense_14/MatMul_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2l
4model_4/sequential_4/dense_12/BiasAdd/ReadVariableOp4model_4/sequential_4/dense_12/BiasAdd/ReadVariableOp2p
6model_4/sequential_4/dense_12/BiasAdd_1/ReadVariableOp6model_4/sequential_4/dense_12/BiasAdd_1/ReadVariableOp2j
3model_4/sequential_4/dense_12/MatMul/ReadVariableOp3model_4/sequential_4/dense_12/MatMul/ReadVariableOp2n
5model_4/sequential_4/dense_12/MatMul_1/ReadVariableOp5model_4/sequential_4/dense_12/MatMul_1/ReadVariableOp2l
4model_4/sequential_4/dense_13/BiasAdd/ReadVariableOp4model_4/sequential_4/dense_13/BiasAdd/ReadVariableOp2p
6model_4/sequential_4/dense_13/BiasAdd_1/ReadVariableOp6model_4/sequential_4/dense_13/BiasAdd_1/ReadVariableOp2j
3model_4/sequential_4/dense_13/MatMul/ReadVariableOp3model_4/sequential_4/dense_13/MatMul/ReadVariableOp2n
5model_4/sequential_4/dense_13/MatMul_1/ReadVariableOp5model_4/sequential_4/dense_13/MatMul_1/ReadVariableOp2l
4model_4/sequential_4/dense_14/BiasAdd/ReadVariableOp4model_4/sequential_4/dense_14/BiasAdd/ReadVariableOp2p
6model_4/sequential_4/dense_14/BiasAdd_1/ReadVariableOp6model_4/sequential_4/dense_14/BiasAdd_1/ReadVariableOp2j
3model_4/sequential_4/dense_14/MatMul/ReadVariableOp3model_4/sequential_4/dense_14/MatMul/ReadVariableOp2n
5model_4/sequential_4/dense_14/MatMul_1/ReadVariableOp5model_4/sequential_4/dense_14/MatMul_1/ReadVariableOp:P L
'
_output_shapes
:         
!
_user_specified_name	input_5
н
L
$__inference__update_step_xla_7019514
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
╔
Ч
*__inference_dense_13_layer_call_fn_7019627

inputs
unknown:22
	unknown_0:2
identityИвStatefulPartitionedCall▀
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
GPU2 *0J 8В *N
fIRG
E__inference_dense_13_layer_call_and_return_conditional_losses_7018952o
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
╗
е
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019129
dense_12_input"
dense_12_7019113:	2
dense_12_7019115:2"
dense_13_7019118:22
dense_13_7019120:2"
dense_14_7019123:2
dense_14_7019125:
identityИв dense_12/StatefulPartitionedCallв dense_13/StatefulPartitionedCallв dense_14/StatefulPartitionedCallА
 dense_12/StatefulPartitionedCallStatefulPartitionedCalldense_12_inputdense_12_7019113dense_12_7019115*
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
GPU2 *0J 8В *N
fIRG
E__inference_dense_12_layer_call_and_return_conditional_losses_7018935Ы
 dense_13/StatefulPartitionedCallStatefulPartitionedCall)dense_12/StatefulPartitionedCall:output:0dense_13_7019118dense_13_7019120*
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
GPU2 *0J 8В *N
fIRG
E__inference_dense_13_layer_call_and_return_conditional_losses_7018952Ы
 dense_14/StatefulPartitionedCallStatefulPartitionedCall)dense_13/StatefulPartitionedCall:output:0dense_14_7019123dense_14_7019125*
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
GPU2 *0J 8В *N
fIRG
E__inference_dense_14_layer_call_and_return_conditional_losses_7018969x
IdentityIdentity)dense_14/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         п
NoOpNoOp!^dense_12/StatefulPartitionedCall!^dense_13/StatefulPartitionedCall!^dense_14/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         	: : : : : : 2D
 dense_12/StatefulPartitionedCall dense_12/StatefulPartitionedCall2D
 dense_13/StatefulPartitionedCall dense_13/StatefulPartitionedCall2D
 dense_14/StatefulPartitionedCall dense_14/StatefulPartitionedCall:W S
'
_output_shapes
:         	
(
_user_specified_namedense_12_input
С	
П
.__inference_sequential_4_layer_call_fn_7019091
dense_12_input
unknown:	2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityИвStatefulPartitionedCallЯ
StatefulPartitionedCallStatefulPartitionedCalldense_12_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
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
GPU2 *0J 8В *R
fMRK
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019059o
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
:         	: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
'
_output_shapes
:         	
(
_user_specified_namedense_12_input
г
Э
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019059

inputs"
dense_12_7019043:	2
dense_12_7019045:2"
dense_13_7019048:22
dense_13_7019050:2"
dense_14_7019053:2
dense_14_7019055:
identityИв dense_12/StatefulPartitionedCallв dense_13/StatefulPartitionedCallв dense_14/StatefulPartitionedCall°
 dense_12/StatefulPartitionedCallStatefulPartitionedCallinputsdense_12_7019043dense_12_7019045*
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
GPU2 *0J 8В *N
fIRG
E__inference_dense_12_layer_call_and_return_conditional_losses_7018935Ы
 dense_13/StatefulPartitionedCallStatefulPartitionedCall)dense_12/StatefulPartitionedCall:output:0dense_13_7019048dense_13_7019050*
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
GPU2 *0J 8В *N
fIRG
E__inference_dense_13_layer_call_and_return_conditional_losses_7018952Ы
 dense_14/StatefulPartitionedCallStatefulPartitionedCall)dense_13/StatefulPartitionedCall:output:0dense_14_7019053dense_14_7019055*
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
GPU2 *0J 8В *N
fIRG
E__inference_dense_14_layer_call_and_return_conditional_losses_7018969x
IdentityIdentity)dense_14/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         п
NoOpNoOp!^dense_12/StatefulPartitionedCall!^dense_13/StatefulPartitionedCall!^dense_14/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         	: : : : : : 2D
 dense_12/StatefulPartitionedCall dense_12/StatefulPartitionedCall2D
 dense_13/StatefulPartitionedCall dense_13/StatefulPartitionedCall2D
 dense_14/StatefulPartitionedCall dense_14/StatefulPartitionedCall:O K
'
_output_shapes
:         	
 
_user_specified_nameinputs
╗
е
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019110
dense_12_input"
dense_12_7019094:	2
dense_12_7019096:2"
dense_13_7019099:22
dense_13_7019101:2"
dense_14_7019104:2
dense_14_7019106:
identityИв dense_12/StatefulPartitionedCallв dense_13/StatefulPartitionedCallв dense_14/StatefulPartitionedCallА
 dense_12/StatefulPartitionedCallStatefulPartitionedCalldense_12_inputdense_12_7019094dense_12_7019096*
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
GPU2 *0J 8В *N
fIRG
E__inference_dense_12_layer_call_and_return_conditional_losses_7018935Ы
 dense_13/StatefulPartitionedCallStatefulPartitionedCall)dense_12/StatefulPartitionedCall:output:0dense_13_7019099dense_13_7019101*
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
GPU2 *0J 8В *N
fIRG
E__inference_dense_13_layer_call_and_return_conditional_losses_7018952Ы
 dense_14/StatefulPartitionedCallStatefulPartitionedCall)dense_13/StatefulPartitionedCall:output:0dense_14_7019104dense_14_7019106*
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
GPU2 *0J 8В *N
fIRG
E__inference_dense_14_layer_call_and_return_conditional_losses_7018969x
IdentityIdentity)dense_14/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         п
NoOpNoOp!^dense_12/StatefulPartitionedCall!^dense_13/StatefulPartitionedCall!^dense_14/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         	: : : : : : 2D
 dense_12/StatefulPartitionedCall dense_12/StatefulPartitionedCall2D
 dense_13/StatefulPartitionedCall dense_13/StatefulPartitionedCall2D
 dense_14/StatefulPartitionedCall dense_14/StatefulPartitionedCall:W S
'
_output_shapes
:         	
(
_user_specified_namedense_12_input
∙
З
.__inference_sequential_4_layer_call_fn_7019548

inputs
unknown:	2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityИвStatefulPartitionedCallЧ
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
GPU2 *0J 8В *R
fMRK
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019059o
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
:         	: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         	
 
_user_specified_nameinputs
щ4
Г

 __inference__traced_save_7019753
file_prefix.
*savev2_dense_12_kernel_read_readvariableop,
(savev2_dense_12_bias_read_readvariableop.
*savev2_dense_13_kernel_read_readvariableop,
(savev2_dense_13_bias_read_readvariableop.
*savev2_dense_14_kernel_read_readvariableop,
(savev2_dense_14_bias_read_readvariableop(
$savev2_iteration_read_readvariableop	,
(savev2_learning_rate_read_readvariableop5
1savev2_adam_m_dense_12_kernel_read_readvariableop5
1savev2_adam_v_dense_12_kernel_read_readvariableop3
/savev2_adam_m_dense_12_bias_read_readvariableop3
/savev2_adam_v_dense_12_bias_read_readvariableop5
1savev2_adam_m_dense_13_kernel_read_readvariableop5
1savev2_adam_v_dense_13_kernel_read_readvariableop3
/savev2_adam_m_dense_13_bias_read_readvariableop3
/savev2_adam_v_dense_13_bias_read_readvariableop5
1savev2_adam_m_dense_14_kernel_read_readvariableop5
1savev2_adam_v_dense_14_kernel_read_readvariableop3
/savev2_adam_m_dense_14_bias_read_readvariableop3
/savev2_adam_v_dense_14_bias_read_readvariableop&
"savev2_total_1_read_readvariableop&
"savev2_count_1_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop
savev2_const

identity_1ИвMergeV2Checkpointsw
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
_temp/partБ
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
value	B : У
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: а

SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*╔	
value┐	B╝	B&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUEB3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPHЯ
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*E
value<B:B B B B B B B B B B B B B B B B B B B B B B B B B ж

SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0*savev2_dense_12_kernel_read_readvariableop(savev2_dense_12_bias_read_readvariableop*savev2_dense_13_kernel_read_readvariableop(savev2_dense_13_bias_read_readvariableop*savev2_dense_14_kernel_read_readvariableop(savev2_dense_14_bias_read_readvariableop$savev2_iteration_read_readvariableop(savev2_learning_rate_read_readvariableop1savev2_adam_m_dense_12_kernel_read_readvariableop1savev2_adam_v_dense_12_kernel_read_readvariableop/savev2_adam_m_dense_12_bias_read_readvariableop/savev2_adam_v_dense_12_bias_read_readvariableop1savev2_adam_m_dense_13_kernel_read_readvariableop1savev2_adam_v_dense_13_kernel_read_readvariableop/savev2_adam_m_dense_13_bias_read_readvariableop/savev2_adam_v_dense_13_bias_read_readvariableop1savev2_adam_m_dense_14_kernel_read_readvariableop1savev2_adam_v_dense_14_kernel_read_readvariableop/savev2_adam_m_dense_14_bias_read_readvariableop/savev2_adam_v_dense_14_bias_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableopsavev2_const"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *'
dtypes
2	Р
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

identity_1Identity_1:output:0*╡
_input_shapesг
а: :	2:2:22:2:2:: : :	2:	2:2:2:22:22:2:2:2:2::: : : : : 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:$ 

_output_shapes

:	2: 
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

:	2:$
 

_output_shapes

:	2: 
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
╣
P
$__inference__update_step_xla_7019509
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
Ь

Ў
E__inference_dense_13_layer_call_and_return_conditional_losses_7019638

inputs0
matmul_readvariableop_resource:22-
biasadd_readvariableop_resource:2
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpt
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
г
Э
I__inference_sequential_4_layer_call_and_return_conditional_losses_7018976

inputs"
dense_12_7018936:	2
dense_12_7018938:2"
dense_13_7018953:22
dense_13_7018955:2"
dense_14_7018970:2
dense_14_7018972:
identityИв dense_12/StatefulPartitionedCallв dense_13/StatefulPartitionedCallв dense_14/StatefulPartitionedCall°
 dense_12/StatefulPartitionedCallStatefulPartitionedCallinputsdense_12_7018936dense_12_7018938*
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
GPU2 *0J 8В *N
fIRG
E__inference_dense_12_layer_call_and_return_conditional_losses_7018935Ы
 dense_13/StatefulPartitionedCallStatefulPartitionedCall)dense_12/StatefulPartitionedCall:output:0dense_13_7018953dense_13_7018955*
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
GPU2 *0J 8В *N
fIRG
E__inference_dense_13_layer_call_and_return_conditional_losses_7018952Ы
 dense_14/StatefulPartitionedCallStatefulPartitionedCall)dense_13/StatefulPartitionedCall:output:0dense_14_7018970dense_14_7018972*
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
GPU2 *0J 8В *N
fIRG
E__inference_dense_14_layer_call_and_return_conditional_losses_7018969x
IdentityIdentity)dense_14/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         п
NoOpNoOp!^dense_12/StatefulPartitionedCall!^dense_13/StatefulPartitionedCall!^dense_14/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         	: : : : : : 2D
 dense_12/StatefulPartitionedCall dense_12/StatefulPartitionedCall2D
 dense_13/StatefulPartitionedCall dense_13/StatefulPartitionedCall2D
 dense_14/StatefulPartitionedCall dense_14/StatefulPartitionedCall:O K
'
_output_shapes
:         	
 
_user_specified_nameinputs
Тf
м
#__inference__traced_restore_7019835
file_prefix2
 assignvariableop_dense_12_kernel:	2.
 assignvariableop_1_dense_12_bias:24
"assignvariableop_2_dense_13_kernel:22.
 assignvariableop_3_dense_13_bias:24
"assignvariableop_4_dense_14_kernel:2.
 assignvariableop_5_dense_14_bias:&
assignvariableop_6_iteration:	 *
 assignvariableop_7_learning_rate: ;
)assignvariableop_8_adam_m_dense_12_kernel:	2;
)assignvariableop_9_adam_v_dense_12_kernel:	26
(assignvariableop_10_adam_m_dense_12_bias:26
(assignvariableop_11_adam_v_dense_12_bias:2<
*assignvariableop_12_adam_m_dense_13_kernel:22<
*assignvariableop_13_adam_v_dense_13_kernel:226
(assignvariableop_14_adam_m_dense_13_bias:26
(assignvariableop_15_adam_v_dense_13_bias:2<
*assignvariableop_16_adam_m_dense_14_kernel:2<
*assignvariableop_17_adam_v_dense_14_kernel:26
(assignvariableop_18_adam_m_dense_14_bias:6
(assignvariableop_19_adam_v_dense_14_bias:%
assignvariableop_20_total_1: %
assignvariableop_21_count_1: #
assignvariableop_22_total: #
assignvariableop_23_count: 
identity_25ИвAssignVariableOpвAssignVariableOp_1вAssignVariableOp_10вAssignVariableOp_11вAssignVariableOp_12вAssignVariableOp_13вAssignVariableOp_14вAssignVariableOp_15вAssignVariableOp_16вAssignVariableOp_17вAssignVariableOp_18вAssignVariableOp_19вAssignVariableOp_2вAssignVariableOp_20вAssignVariableOp_21вAssignVariableOp_22вAssignVariableOp_23вAssignVariableOp_3вAssignVariableOp_4вAssignVariableOp_5вAssignVariableOp_6вAssignVariableOp_7вAssignVariableOp_8вAssignVariableOp_9г

RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*╔	
value┐	B╝	B&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUEB3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPHв
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*E
value<B:B B B B B B B B B B B B B B B B B B B B B B B B B Ы
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*x
_output_shapesf
d:::::::::::::::::::::::::*'
dtypes
2	[
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:│
AssignVariableOpAssignVariableOp assignvariableop_dense_12_kernelIdentity:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:╖
AssignVariableOp_1AssignVariableOp assignvariableop_1_dense_12_biasIdentity_1:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:╣
AssignVariableOp_2AssignVariableOp"assignvariableop_2_dense_13_kernelIdentity_2:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:╖
AssignVariableOp_3AssignVariableOp assignvariableop_3_dense_13_biasIdentity_3:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:╣
AssignVariableOp_4AssignVariableOp"assignvariableop_4_dense_14_kernelIdentity_4:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:╖
AssignVariableOp_5AssignVariableOp assignvariableop_5_dense_14_biasIdentity_5:output:0"/device:CPU:0*&
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
:╖
AssignVariableOp_7AssignVariableOp assignvariableop_7_learning_rateIdentity_7:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:└
AssignVariableOp_8AssignVariableOp)assignvariableop_8_adam_m_dense_12_kernelIdentity_8:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:└
AssignVariableOp_9AssignVariableOp)assignvariableop_9_adam_v_dense_12_kernelIdentity_9:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:┴
AssignVariableOp_10AssignVariableOp(assignvariableop_10_adam_m_dense_12_biasIdentity_10:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:┴
AssignVariableOp_11AssignVariableOp(assignvariableop_11_adam_v_dense_12_biasIdentity_11:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:├
AssignVariableOp_12AssignVariableOp*assignvariableop_12_adam_m_dense_13_kernelIdentity_12:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:├
AssignVariableOp_13AssignVariableOp*assignvariableop_13_adam_v_dense_13_kernelIdentity_13:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:┴
AssignVariableOp_14AssignVariableOp(assignvariableop_14_adam_m_dense_13_biasIdentity_14:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:┴
AssignVariableOp_15AssignVariableOp(assignvariableop_15_adam_v_dense_13_biasIdentity_15:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:├
AssignVariableOp_16AssignVariableOp*assignvariableop_16_adam_m_dense_14_kernelIdentity_16:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:├
AssignVariableOp_17AssignVariableOp*assignvariableop_17_adam_v_dense_14_kernelIdentity_17:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:┴
AssignVariableOp_18AssignVariableOp(assignvariableop_18_adam_m_dense_14_biasIdentity_18:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:┴
AssignVariableOp_19AssignVariableOp(assignvariableop_19_adam_v_dense_14_biasIdentity_19:output:0"/device:CPU:0*&
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
_user_specified_namefile_prefix
╔
Ч
*__inference_dense_14_layer_call_fn_7019647

inputs
unknown:2
	unknown_0:
identityИвStatefulPartitionedCall▀
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
GPU2 *0J 8В *N
fIRG
E__inference_dense_14_layer_call_and_return_conditional_losses_7018969o
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
┐
Ш
D__inference_model_4_layer_call_and_return_conditional_losses_7019298
input_5&
sequential_4_7019276:	2"
sequential_4_7019278:2&
sequential_4_7019280:22"
sequential_4_7019282:2&
sequential_4_7019284:2"
sequential_4_7019286:
identityИв$sequential_4/StatefulPartitionedCallв&sequential_4/StatefulPartitionedCall_1
.tf.__operators__.getitem_8/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"    	   Б
0tf.__operators__.getitem_8/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        Б
0tf.__operators__.getitem_8/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ╤
(tf.__operators__.getitem_8/strided_sliceStridedSliceinput_57tf.__operators__.getitem_8/strided_slice/stack:output:09tf.__operators__.getitem_8/strided_slice/stack_1:output:09tf.__operators__.getitem_8/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         	*

begin_mask*
end_mask
.tf.__operators__.getitem_7/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        Б
0tf.__operators__.getitem_7/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"    	   Б
0tf.__operators__.getitem_7/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ╤
(tf.__operators__.getitem_7/strided_sliceStridedSliceinput_57tf.__operators__.getitem_7/strided_slice/stack:output:09tf.__operators__.getitem_7/strided_slice/stack_1:output:09tf.__operators__.getitem_7/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         	*

begin_mask*
end_maskУ
$sequential_4/StatefulPartitionedCallStatefulPartitionedCall1tf.__operators__.getitem_7/strided_slice:output:0sequential_4_7019276sequential_4_7019278sequential_4_7019280sequential_4_7019282sequential_4_7019284sequential_4_7019286*
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
GPU2 *0J 8В *R
fMRK
I__inference_sequential_4_layer_call_and_return_conditional_losses_7018976Х
&sequential_4/StatefulPartitionedCall_1StatefulPartitionedCall1tf.__operators__.getitem_8/strided_slice:output:0sequential_4_7019276sequential_4_7019278sequential_4_7019280sequential_4_7019282sequential_4_7019284sequential_4_7019286*
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
GPU2 *0J 8В *R
fMRK
I__inference_sequential_4_layer_call_and_return_conditional_losses_7018976├
tf.stack_3/stackPack-sequential_4/StatefulPartitionedCall:output:0/sequential_4/StatefulPartitionedCall_1:output:0*
N*
T0*+
_output_shapes
:         *

axisl
IdentityIdentitytf.stack_3/stack:output:0^NoOp*
T0*+
_output_shapes
:         Ц
NoOpNoOp%^sequential_4/StatefulPartitionedCall'^sequential_4/StatefulPartitionedCall_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2L
$sequential_4/StatefulPartitionedCall$sequential_4/StatefulPartitionedCall2P
&sequential_4/StatefulPartitionedCall_1&sequential_4/StatefulPartitionedCall_1:P L
'
_output_shapes
:         
!
_user_specified_name	input_5
╔
Ч
*__inference_dense_12_layer_call_fn_7019607

inputs
unknown:	2
	unknown_0:2
identityИвStatefulPartitionedCall▀
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
GPU2 *0J 8В *N
fIRG
E__inference_dense_12_layer_call_and_return_conditional_losses_7018935o
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
:         	: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         	
 
_user_specified_nameinputs
е
К
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019573

inputs9
'dense_12_matmul_readvariableop_resource:	26
(dense_12_biasadd_readvariableop_resource:29
'dense_13_matmul_readvariableop_resource:226
(dense_13_biasadd_readvariableop_resource:29
'dense_14_matmul_readvariableop_resource:26
(dense_14_biasadd_readvariableop_resource:
identityИвdense_12/BiasAdd/ReadVariableOpвdense_12/MatMul/ReadVariableOpвdense_13/BiasAdd/ReadVariableOpвdense_13/MatMul/ReadVariableOpвdense_14/BiasAdd/ReadVariableOpвdense_14/MatMul/ReadVariableOpЖ
dense_12/MatMul/ReadVariableOpReadVariableOp'dense_12_matmul_readvariableop_resource*
_output_shapes

:	2*
dtype0{
dense_12/MatMulMatMulinputs&dense_12/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2Д
dense_12/BiasAdd/ReadVariableOpReadVariableOp(dense_12_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0С
dense_12/BiasAddBiasAdddense_12/MatMul:product:0'dense_12/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2b
dense_12/ReluReludense_12/BiasAdd:output:0*
T0*'
_output_shapes
:         2Ж
dense_13/MatMul/ReadVariableOpReadVariableOp'dense_13_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Р
dense_13/MatMulMatMuldense_12/Relu:activations:0&dense_13/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2Д
dense_13/BiasAdd/ReadVariableOpReadVariableOp(dense_13_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0С
dense_13/BiasAddBiasAdddense_13/MatMul:product:0'dense_13/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2b
dense_13/ReluReludense_13/BiasAdd:output:0*
T0*'
_output_shapes
:         2Ж
dense_14/MatMul/ReadVariableOpReadVariableOp'dense_14_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0Р
dense_14/MatMulMatMuldense_13/Relu:activations:0&dense_14/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         Д
dense_14/BiasAdd/ReadVariableOpReadVariableOp(dense_14_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0С
dense_14/BiasAddBiasAdddense_14/MatMul:product:0'dense_14/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         h
dense_14/SigmoidSigmoiddense_14/BiasAdd:output:0*
T0*'
_output_shapes
:         c
IdentityIdentitydense_14/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:         П
NoOpNoOp ^dense_12/BiasAdd/ReadVariableOp^dense_12/MatMul/ReadVariableOp ^dense_13/BiasAdd/ReadVariableOp^dense_13/MatMul/ReadVariableOp ^dense_14/BiasAdd/ReadVariableOp^dense_14/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         	: : : : : : 2B
dense_12/BiasAdd/ReadVariableOpdense_12/BiasAdd/ReadVariableOp2@
dense_12/MatMul/ReadVariableOpdense_12/MatMul/ReadVariableOp2B
dense_13/BiasAdd/ReadVariableOpdense_13/BiasAdd/ReadVariableOp2@
dense_13/MatMul/ReadVariableOpdense_13/MatMul/ReadVariableOp2B
dense_14/BiasAdd/ReadVariableOpdense_14/BiasAdd/ReadVariableOp2@
dense_14/MatMul/ReadVariableOpdense_14/MatMul/ReadVariableOp:O K
'
_output_shapes
:         	
 
_user_specified_nameinputs
Ы

Ў
E__inference_dense_14_layer_call_and_return_conditional_losses_7019658

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpt
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
:         V
SigmoidSigmoidBiasAdd:output:0*
T0*'
_output_shapes
:         Z
IdentityIdentitySigmoid:y:0^NoOp*
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
Ы

Ў
E__inference_dense_14_layer_call_and_return_conditional_losses_7018969

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpt
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
:         V
SigmoidSigmoidBiasAdd:output:0*
T0*'
_output_shapes
:         Z
IdentityIdentitySigmoid:y:0^NoOp*
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
╣
P
$__inference__update_step_xla_7019499
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
ЖM
─
D__inference_model_4_layer_call_and_return_conditional_losses_7019435

inputsF
4sequential_4_dense_12_matmul_readvariableop_resource:	2C
5sequential_4_dense_12_biasadd_readvariableop_resource:2F
4sequential_4_dense_13_matmul_readvariableop_resource:22C
5sequential_4_dense_13_biasadd_readvariableop_resource:2F
4sequential_4_dense_14_matmul_readvariableop_resource:2C
5sequential_4_dense_14_biasadd_readvariableop_resource:
identityИв,sequential_4/dense_12/BiasAdd/ReadVariableOpв.sequential_4/dense_12/BiasAdd_1/ReadVariableOpв+sequential_4/dense_12/MatMul/ReadVariableOpв-sequential_4/dense_12/MatMul_1/ReadVariableOpв,sequential_4/dense_13/BiasAdd/ReadVariableOpв.sequential_4/dense_13/BiasAdd_1/ReadVariableOpв+sequential_4/dense_13/MatMul/ReadVariableOpв-sequential_4/dense_13/MatMul_1/ReadVariableOpв,sequential_4/dense_14/BiasAdd/ReadVariableOpв.sequential_4/dense_14/BiasAdd_1/ReadVariableOpв+sequential_4/dense_14/MatMul/ReadVariableOpв-sequential_4/dense_14/MatMul_1/ReadVariableOp
.tf.__operators__.getitem_8/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"    	   Б
0tf.__operators__.getitem_8/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        Б
0tf.__operators__.getitem_8/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ╨
(tf.__operators__.getitem_8/strided_sliceStridedSliceinputs7tf.__operators__.getitem_8/strided_slice/stack:output:09tf.__operators__.getitem_8/strided_slice/stack_1:output:09tf.__operators__.getitem_8/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         	*

begin_mask*
end_mask
.tf.__operators__.getitem_7/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        Б
0tf.__operators__.getitem_7/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"    	   Б
0tf.__operators__.getitem_7/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ╨
(tf.__operators__.getitem_7/strided_sliceStridedSliceinputs7tf.__operators__.getitem_7/strided_slice/stack:output:09tf.__operators__.getitem_7/strided_slice/stack_1:output:09tf.__operators__.getitem_7/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         	*

begin_mask*
end_maskа
+sequential_4/dense_12/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_12_matmul_readvariableop_resource*
_output_shapes

:	2*
dtype0└
sequential_4/dense_12/MatMulMatMul1tf.__operators__.getitem_7/strided_slice:output:03sequential_4/dense_12/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2Ю
,sequential_4/dense_12/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_12_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0╕
sequential_4/dense_12/BiasAddBiasAdd&sequential_4/dense_12/MatMul:product:04sequential_4/dense_12/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2|
sequential_4/dense_12/ReluRelu&sequential_4/dense_12/BiasAdd:output:0*
T0*'
_output_shapes
:         2а
+sequential_4/dense_13/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_13_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0╖
sequential_4/dense_13/MatMulMatMul(sequential_4/dense_12/Relu:activations:03sequential_4/dense_13/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2Ю
,sequential_4/dense_13/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_13_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0╕
sequential_4/dense_13/BiasAddBiasAdd&sequential_4/dense_13/MatMul:product:04sequential_4/dense_13/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2|
sequential_4/dense_13/ReluRelu&sequential_4/dense_13/BiasAdd:output:0*
T0*'
_output_shapes
:         2а
+sequential_4/dense_14/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_14_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0╖
sequential_4/dense_14/MatMulMatMul(sequential_4/dense_13/Relu:activations:03sequential_4/dense_14/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         Ю
,sequential_4/dense_14/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_14_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0╕
sequential_4/dense_14/BiasAddBiasAdd&sequential_4/dense_14/MatMul:product:04sequential_4/dense_14/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         В
sequential_4/dense_14/SigmoidSigmoid&sequential_4/dense_14/BiasAdd:output:0*
T0*'
_output_shapes
:         в
-sequential_4/dense_12/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_12_matmul_readvariableop_resource*
_output_shapes

:	2*
dtype0─
sequential_4/dense_12/MatMul_1MatMul1tf.__operators__.getitem_8/strided_slice:output:05sequential_4/dense_12/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2а
.sequential_4/dense_12/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_12_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0╛
sequential_4/dense_12/BiasAdd_1BiasAdd(sequential_4/dense_12/MatMul_1:product:06sequential_4/dense_12/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2А
sequential_4/dense_12/Relu_1Relu(sequential_4/dense_12/BiasAdd_1:output:0*
T0*'
_output_shapes
:         2в
-sequential_4/dense_13/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_13_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0╜
sequential_4/dense_13/MatMul_1MatMul*sequential_4/dense_12/Relu_1:activations:05sequential_4/dense_13/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2а
.sequential_4/dense_13/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_13_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0╛
sequential_4/dense_13/BiasAdd_1BiasAdd(sequential_4/dense_13/MatMul_1:product:06sequential_4/dense_13/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2А
sequential_4/dense_13/Relu_1Relu(sequential_4/dense_13/BiasAdd_1:output:0*
T0*'
_output_shapes
:         2в
-sequential_4/dense_14/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_14_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0╜
sequential_4/dense_14/MatMul_1MatMul*sequential_4/dense_13/Relu_1:activations:05sequential_4/dense_14/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         а
.sequential_4/dense_14/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_14_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0╛
sequential_4/dense_14/BiasAdd_1BiasAdd(sequential_4/dense_14/MatMul_1:product:06sequential_4/dense_14/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:         Ж
sequential_4/dense_14/Sigmoid_1Sigmoid(sequential_4/dense_14/BiasAdd_1:output:0*
T0*'
_output_shapes
:         л
tf.stack_3/stackPack!sequential_4/dense_14/Sigmoid:y:0#sequential_4/dense_14/Sigmoid_1:y:0*
N*
T0*+
_output_shapes
:         *

axisl
IdentityIdentitytf.stack_3/stack:output:0^NoOp*
T0*+
_output_shapes
:         А
NoOpNoOp-^sequential_4/dense_12/BiasAdd/ReadVariableOp/^sequential_4/dense_12/BiasAdd_1/ReadVariableOp,^sequential_4/dense_12/MatMul/ReadVariableOp.^sequential_4/dense_12/MatMul_1/ReadVariableOp-^sequential_4/dense_13/BiasAdd/ReadVariableOp/^sequential_4/dense_13/BiasAdd_1/ReadVariableOp,^sequential_4/dense_13/MatMul/ReadVariableOp.^sequential_4/dense_13/MatMul_1/ReadVariableOp-^sequential_4/dense_14/BiasAdd/ReadVariableOp/^sequential_4/dense_14/BiasAdd_1/ReadVariableOp,^sequential_4/dense_14/MatMul/ReadVariableOp.^sequential_4/dense_14/MatMul_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2\
,sequential_4/dense_12/BiasAdd/ReadVariableOp,sequential_4/dense_12/BiasAdd/ReadVariableOp2`
.sequential_4/dense_12/BiasAdd_1/ReadVariableOp.sequential_4/dense_12/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_12/MatMul/ReadVariableOp+sequential_4/dense_12/MatMul/ReadVariableOp2^
-sequential_4/dense_12/MatMul_1/ReadVariableOp-sequential_4/dense_12/MatMul_1/ReadVariableOp2\
,sequential_4/dense_13/BiasAdd/ReadVariableOp,sequential_4/dense_13/BiasAdd/ReadVariableOp2`
.sequential_4/dense_13/BiasAdd_1/ReadVariableOp.sequential_4/dense_13/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_13/MatMul/ReadVariableOp+sequential_4/dense_13/MatMul/ReadVariableOp2^
-sequential_4/dense_13/MatMul_1/ReadVariableOp-sequential_4/dense_13/MatMul_1/ReadVariableOp2\
,sequential_4/dense_14/BiasAdd/ReadVariableOp,sequential_4/dense_14/BiasAdd/ReadVariableOp2`
.sequential_4/dense_14/BiasAdd_1/ReadVariableOp.sequential_4/dense_14/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_14/MatMul/ReadVariableOp+sequential_4/dense_14/MatMul/ReadVariableOp2^
-sequential_4/dense_14/MatMul_1/ReadVariableOp-sequential_4/dense_14/MatMul_1/ReadVariableOp:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
∙
З
.__inference_sequential_4_layer_call_fn_7019531

inputs
unknown:	2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityИвStatefulPartitionedCallЧ
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
GPU2 *0J 8В *R
fMRK
I__inference_sequential_4_layer_call_and_return_conditional_losses_7018976o
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
:         	: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         	
 
_user_specified_nameinputs
ў
В
)__inference_model_4_layer_call_fn_7019369

inputs
unknown:	2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
identityИвStatefulPartitionedCallЦ
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
GPU2 *0J 8В *M
fHRF
D__inference_model_4_layer_call_and_return_conditional_losses_7019166s
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
:         : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
╗
Ч
D__inference_model_4_layer_call_and_return_conditional_losses_7019166

inputs&
sequential_4_7019144:	2"
sequential_4_7019146:2&
sequential_4_7019148:22"
sequential_4_7019150:2&
sequential_4_7019152:2"
sequential_4_7019154:
identityИв$sequential_4/StatefulPartitionedCallв&sequential_4/StatefulPartitionedCall_1
.tf.__operators__.getitem_8/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"    	   Б
0tf.__operators__.getitem_8/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        Б
0tf.__operators__.getitem_8/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ╨
(tf.__operators__.getitem_8/strided_sliceStridedSliceinputs7tf.__operators__.getitem_8/strided_slice/stack:output:09tf.__operators__.getitem_8/strided_slice/stack_1:output:09tf.__operators__.getitem_8/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         	*

begin_mask*
end_mask
.tf.__operators__.getitem_7/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        Б
0tf.__operators__.getitem_7/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"    	   Б
0tf.__operators__.getitem_7/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ╨
(tf.__operators__.getitem_7/strided_sliceStridedSliceinputs7tf.__operators__.getitem_7/strided_slice/stack:output:09tf.__operators__.getitem_7/strided_slice/stack_1:output:09tf.__operators__.getitem_7/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         	*

begin_mask*
end_maskУ
$sequential_4/StatefulPartitionedCallStatefulPartitionedCall1tf.__operators__.getitem_7/strided_slice:output:0sequential_4_7019144sequential_4_7019146sequential_4_7019148sequential_4_7019150sequential_4_7019152sequential_4_7019154*
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
GPU2 *0J 8В *R
fMRK
I__inference_sequential_4_layer_call_and_return_conditional_losses_7018976Х
&sequential_4/StatefulPartitionedCall_1StatefulPartitionedCall1tf.__operators__.getitem_8/strided_slice:output:0sequential_4_7019144sequential_4_7019146sequential_4_7019148sequential_4_7019150sequential_4_7019152sequential_4_7019154*
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
GPU2 *0J 8В *R
fMRK
I__inference_sequential_4_layer_call_and_return_conditional_losses_7018976├
tf.stack_3/stackPack-sequential_4/StatefulPartitionedCall:output:0/sequential_4/StatefulPartitionedCall_1:output:0*
N*
T0*+
_output_shapes
:         *

axisl
IdentityIdentitytf.stack_3/stack:output:0^NoOp*
T0*+
_output_shapes
:         Ц
NoOpNoOp%^sequential_4/StatefulPartitionedCall'^sequential_4/StatefulPartitionedCall_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2L
$sequential_4/StatefulPartitionedCall$sequential_4/StatefulPartitionedCall2P
&sequential_4/StatefulPartitionedCall_1&sequential_4/StatefulPartitionedCall_1:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
╗
Ч
D__inference_model_4_layer_call_and_return_conditional_losses_7019233

inputs&
sequential_4_7019211:	2"
sequential_4_7019213:2&
sequential_4_7019215:22"
sequential_4_7019217:2&
sequential_4_7019219:2"
sequential_4_7019221:
identityИв$sequential_4/StatefulPartitionedCallв&sequential_4/StatefulPartitionedCall_1
.tf.__operators__.getitem_8/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"    	   Б
0tf.__operators__.getitem_8/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        Б
0tf.__operators__.getitem_8/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ╨
(tf.__operators__.getitem_8/strided_sliceStridedSliceinputs7tf.__operators__.getitem_8/strided_slice/stack:output:09tf.__operators__.getitem_8/strided_slice/stack_1:output:09tf.__operators__.getitem_8/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         	*

begin_mask*
end_mask
.tf.__operators__.getitem_7/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        Б
0tf.__operators__.getitem_7/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"    	   Б
0tf.__operators__.getitem_7/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      ╨
(tf.__operators__.getitem_7/strided_sliceStridedSliceinputs7tf.__operators__.getitem_7/strided_slice/stack:output:09tf.__operators__.getitem_7/strided_slice/stack_1:output:09tf.__operators__.getitem_7/strided_slice/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         	*

begin_mask*
end_maskУ
$sequential_4/StatefulPartitionedCallStatefulPartitionedCall1tf.__operators__.getitem_7/strided_slice:output:0sequential_4_7019211sequential_4_7019213sequential_4_7019215sequential_4_7019217sequential_4_7019219sequential_4_7019221*
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
GPU2 *0J 8В *R
fMRK
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019059Х
&sequential_4/StatefulPartitionedCall_1StatefulPartitionedCall1tf.__operators__.getitem_8/strided_slice:output:0sequential_4_7019211sequential_4_7019213sequential_4_7019215sequential_4_7019217sequential_4_7019219sequential_4_7019221*
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
GPU2 *0J 8В *R
fMRK
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019059├
tf.stack_3/stackPack-sequential_4/StatefulPartitionedCall:output:0/sequential_4/StatefulPartitionedCall_1:output:0*
N*
T0*+
_output_shapes
:         *

axisl
IdentityIdentitytf.stack_3/stack:output:0^NoOp*
T0*+
_output_shapes
:         Ц
NoOpNoOp%^sequential_4/StatefulPartitionedCall'^sequential_4/StatefulPartitionedCall_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : : : : : 2L
$sequential_4/StatefulPartitionedCall$sequential_4/StatefulPartitionedCall2P
&sequential_4/StatefulPartitionedCall_1&sequential_4/StatefulPartitionedCall_1:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
н
L
$__inference__update_step_xla_7019494
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
variable"Ж
L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*▒
serving_defaultЭ
;
input_50
serving_default_input_5:0         B

tf.stack_34
StatefulPartitionedCall:0         tensorflow/serving/predict:▄╡
╛
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
Я
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
┘
&trace_0
'trace_1
(trace_2
)trace_32ю
)__inference_model_4_layer_call_fn_7019181
)__inference_model_4_layer_call_fn_7019369
)__inference_model_4_layer_call_fn_7019386
)__inference_model_4_layer_call_fn_7019265┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 z&trace_0z'trace_1z(trace_2z)trace_3
┼
*trace_0
+trace_1
,trace_2
-trace_32┌
D__inference_model_4_layer_call_and_return_conditional_losses_7019435
D__inference_model_4_layer_call_and_return_conditional_losses_7019484
D__inference_model_4_layer_call_and_return_conditional_losses_7019298
D__inference_model_4_layer_call_and_return_conditional_losses_7019331┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 z*trace_0z+trace_1z,trace_2z-trace_3
═B╩
"__inference__wrapped_model_7018917input_5"Ш
С▓Н
FullArgSpec
argsЪ 
varargsjargs
varkwjkwargs
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
Ь
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
н
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
Ptrace_32В
.__inference_sequential_4_layer_call_fn_7018991
.__inference_sequential_4_layer_call_fn_7019531
.__inference_sequential_4_layer_call_fn_7019548
.__inference_sequential_4_layer_call_fn_7019091┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 zMtrace_0zNtrace_1zOtrace_2zPtrace_3
┘
Qtrace_0
Rtrace_1
Strace_2
Ttrace_32ю
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019573
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019598
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019110
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019129┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 zQtrace_0zRtrace_1zStrace_2zTtrace_3
"
_generic_user_object
!:	22dense_12/kernel
:22dense_12/bias
!:222dense_13/kernel
:22dense_13/bias
!:22dense_14/kernel
:2dense_14/bias
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
√B°
)__inference_model_4_layer_call_fn_7019181input_5"┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
·Bў
)__inference_model_4_layer_call_fn_7019369inputs"┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
·Bў
)__inference_model_4_layer_call_fn_7019386inputs"┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
√B°
)__inference_model_4_layer_call_fn_7019265input_5"┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ХBТ
D__inference_model_4_layer_call_and_return_conditional_losses_7019435inputs"┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ХBТ
D__inference_model_4_layer_call_and_return_conditional_losses_7019484inputs"┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ЦBУ
D__inference_model_4_layer_call_and_return_conditional_losses_7019298input_5"┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ЦBУ
D__inference_model_4_layer_call_and_return_conditional_losses_7019331input_5"┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
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
┐
ctrace_0
dtrace_1
etrace_2
ftrace_3
gtrace_4
htrace_52а
$__inference__update_step_xla_7019489
$__inference__update_step_xla_7019494
$__inference__update_step_xla_7019499
$__inference__update_step_xla_7019504
$__inference__update_step_xla_7019509
$__inference__update_step_xla_7019514╣
о▓к
FullArgSpec2
args*Ъ'
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

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 0zctrace_0zdtrace_1zetrace_2zftrace_3zgtrace_4zhtrace_5
╠B╔
%__inference_signature_wrapper_7019352input_5"Ф
Н▓Й
FullArgSpec
argsЪ 
varargs
 
varkwjkwargs
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
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
н
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
ю
ntrace_02╤
*__inference_dense_12_layer_call_fn_7019607в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 zntrace_0
Й
otrace_02ь
E__inference_dense_12_layer_call_and_return_conditional_losses_7019618в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
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
н
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
ю
utrace_02╤
*__inference_dense_13_layer_call_fn_7019627в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 zutrace_0
Й
vtrace_02ь
E__inference_dense_13_layer_call_and_return_conditional_losses_7019638в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
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
н
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
ю
|trace_02╤
*__inference_dense_14_layer_call_fn_7019647в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 z|trace_0
Й
}trace_02ь
E__inference_dense_14_layer_call_and_return_conditional_losses_7019658в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
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
ЗBД
.__inference_sequential_4_layer_call_fn_7018991dense_12_input"┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
 B№
.__inference_sequential_4_layer_call_fn_7019531inputs"┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
 B№
.__inference_sequential_4_layer_call_fn_7019548inputs"┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ЗBД
.__inference_sequential_4_layer_call_fn_7019091dense_12_input"┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ЪBЧ
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019573inputs"┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ЪBЧ
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019598inputs"┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
вBЯ
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019110dense_12_input"┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
вBЯ
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019129dense_12_input"┐
╢▓▓
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
P
~	variables
	keras_api

Аtotal

Бcount"
_tf_keras_metric
c
В	variables
Г	keras_api

Дtotal

Еcount
Ж
_fn_kwargs"
_tf_keras_metric
&:$	22Adam/m/dense_12/kernel
&:$	22Adam/v/dense_12/kernel
 :22Adam/m/dense_12/bias
 :22Adam/v/dense_12/bias
&:$222Adam/m/dense_13/kernel
&:$222Adam/v/dense_13/kernel
 :22Adam/m/dense_13/bias
 :22Adam/v/dense_13/bias
&:$22Adam/m/dense_14/kernel
&:$22Adam/v/dense_14/kernel
 :2Adam/m/dense_14/bias
 :2Adam/v/dense_14/bias
∙BЎ
$__inference__update_step_xla_7019489gradientvariable"╖
о▓к
FullArgSpec2
args*Ъ'
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

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
∙BЎ
$__inference__update_step_xla_7019494gradientvariable"╖
о▓к
FullArgSpec2
args*Ъ'
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

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
∙BЎ
$__inference__update_step_xla_7019499gradientvariable"╖
о▓к
FullArgSpec2
args*Ъ'
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

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
∙BЎ
$__inference__update_step_xla_7019504gradientvariable"╖
о▓к
FullArgSpec2
args*Ъ'
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

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
∙BЎ
$__inference__update_step_xla_7019509gradientvariable"╖
о▓к
FullArgSpec2
args*Ъ'
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

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
∙BЎ
$__inference__update_step_xla_7019514gradientvariable"╖
о▓к
FullArgSpec2
args*Ъ'
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

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
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
▐B█
*__inference_dense_12_layer_call_fn_7019607inputs"в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
∙BЎ
E__inference_dense_12_layer_call_and_return_conditional_losses_7019618inputs"в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
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
▐B█
*__inference_dense_13_layer_call_fn_7019627inputs"в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
∙BЎ
E__inference_dense_13_layer_call_and_return_conditional_losses_7019638inputs"в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
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
▐B█
*__inference_dense_14_layer_call_fn_7019647inputs"в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
∙BЎ
E__inference_dense_14_layer_call_and_return_conditional_losses_7019658inputs"в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
0
А0
Б1"
trackable_list_wrapper
-
~	variables"
_generic_user_object
:  (2total
:  (2count
0
Д0
Е1"
trackable_list_wrapper
.
В	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapperЦ
$__inference__update_step_xla_7019489nhвe
^в[
К
gradient	2
4Т1	в
·	2
А
p
` VariableSpec 
`└ЛБ─╞я?
к "
 О
$__inference__update_step_xla_7019494f`в]
VвS
К
gradient2
0Т-	в
·2
А
p
` VariableSpec 
`А╚А─╞я?
к "
 Ц
$__inference__update_step_xla_7019499nhвe
^в[
К
gradient22
4Т1	в
·22
А
p
` VariableSpec 
`└о╧─╞я?
к "
 О
$__inference__update_step_xla_7019504f`в]
VвS
К
gradient2
0Т-	в
·2
А
p
` VariableSpec 
`аП╧─╞я?
к "
 Ц
$__inference__update_step_xla_7019509nhвe
^в[
К
gradient2
4Т1	в
·2
А
p
` VariableSpec 
`ач╬─╞я?
к "
 О
$__inference__update_step_xla_7019514f`в]
VвS
К
gradient
0Т-	в
·
А
p
` VariableSpec 
`└Є╬─╞я?
к "
 Э
"__inference__wrapped_model_7018917w 0в-
&в#
!К
input_5         
к ";к8
6

tf.stack_3(К%

tf_stack_3         м
E__inference_dense_12_layer_call_and_return_conditional_losses_7019618c/в,
%в"
 К
inputs         	
к ",в)
"К
tensor_0         2
Ъ Ж
*__inference_dense_12_layer_call_fn_7019607X/в,
%в"
 К
inputs         	
к "!К
unknown         2м
E__inference_dense_13_layer_call_and_return_conditional_losses_7019638c/в,
%в"
 К
inputs         2
к ",в)
"К
tensor_0         2
Ъ Ж
*__inference_dense_13_layer_call_fn_7019627X/в,
%в"
 К
inputs         2
к "!К
unknown         2м
E__inference_dense_14_layer_call_and_return_conditional_losses_7019658c /в,
%в"
 К
inputs         2
к ",в)
"К
tensor_0         
Ъ Ж
*__inference_dense_14_layer_call_fn_7019647X /в,
%в"
 К
inputs         2
к "!К
unknown         ╝
D__inference_model_4_layer_call_and_return_conditional_losses_7019298t 8в5
.в+
!К
input_5         
p 

 
к "0в-
&К#
tensor_0         
Ъ ╝
D__inference_model_4_layer_call_and_return_conditional_losses_7019331t 8в5
.в+
!К
input_5         
p

 
к "0в-
&К#
tensor_0         
Ъ ╗
D__inference_model_4_layer_call_and_return_conditional_losses_7019435s 7в4
-в*
 К
inputs         
p 

 
к "0в-
&К#
tensor_0         
Ъ ╗
D__inference_model_4_layer_call_and_return_conditional_losses_7019484s 7в4
-в*
 К
inputs         
p

 
к "0в-
&К#
tensor_0         
Ъ Ц
)__inference_model_4_layer_call_fn_7019181i 8в5
.в+
!К
input_5         
p 

 
к "%К"
unknown         Ц
)__inference_model_4_layer_call_fn_7019265i 8в5
.в+
!К
input_5         
p

 
к "%К"
unknown         Х
)__inference_model_4_layer_call_fn_7019369h 7в4
-в*
 К
inputs         
p 

 
к "%К"
unknown         Х
)__inference_model_4_layer_call_fn_7019386h 7в4
-в*
 К
inputs         
p

 
к "%К"
unknown         ─
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019110w ?в<
5в2
(К%
dense_12_input         	
p 

 
к ",в)
"К
tensor_0         
Ъ ─
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019129w ?в<
5в2
(К%
dense_12_input         	
p

 
к ",в)
"К
tensor_0         
Ъ ╝
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019573o 7в4
-в*
 К
inputs         	
p 

 
к ",в)
"К
tensor_0         
Ъ ╝
I__inference_sequential_4_layer_call_and_return_conditional_losses_7019598o 7в4
-в*
 К
inputs         	
p

 
к ",в)
"К
tensor_0         
Ъ Ю
.__inference_sequential_4_layer_call_fn_7018991l ?в<
5в2
(К%
dense_12_input         	
p 

 
к "!К
unknown         Ю
.__inference_sequential_4_layer_call_fn_7019091l ?в<
5в2
(К%
dense_12_input         	
p

 
к "!К
unknown         Ц
.__inference_sequential_4_layer_call_fn_7019531d 7в4
-в*
 К
inputs         	
p 

 
к "!К
unknown         Ц
.__inference_sequential_4_layer_call_fn_7019548d 7в4
-в*
 К
inputs         	
p

 
к "!К
unknown         м
%__inference_signature_wrapper_7019352В ;в8
в 
1к.
,
input_5!К
input_5         ";к8
6

tf.stack_3(К%

tf_stack_3         