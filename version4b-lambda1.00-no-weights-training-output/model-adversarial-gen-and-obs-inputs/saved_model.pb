шщ
ЄЫ
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
h
ConcatV2
values"T*N
axis"Tidx
output"T"
Nint(0"	
Ttype"
Tidxtype0:
2	
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
9
	IdentityN

input2T
output2T"
T
list(type)(0
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
Ѕ
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
executor_typestring И®
@
StaticRegexFullMatch	
input

output
"
patternstring
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
 И"serve*2.11.02unknown8ОО
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
b
count_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_2
[
count_2/Read/ReadVariableOpReadVariableOpcount_2*
_output_shapes
: *
dtype0
b
total_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_2
[
total_2/Read/ReadVariableOpReadVariableOptotal_2*
_output_shapes
: *
dtype0
v
Adam/v/Adv/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_nameAdam/v/Adv/bias
o
#Adam/v/Adv/bias/Read/ReadVariableOpReadVariableOpAdam/v/Adv/bias*
_output_shapes
:*
dtype0
v
Adam/m/Adv/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_nameAdam/m/Adv/bias
o
#Adam/m/Adv/bias/Read/ReadVariableOpReadVariableOpAdam/m/Adv/bias*
_output_shapes
:*
dtype0
~
Adam/v/Adv/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*"
shared_nameAdam/v/Adv/kernel
w
%Adam/v/Adv/kernel/Read/ReadVariableOpReadVariableOpAdam/v/Adv/kernel*
_output_shapes

:2*
dtype0
~
Adam/m/Adv/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*"
shared_nameAdam/m/Adv/kernel
w
%Adam/m/Adv/kernel/Read/ReadVariableOpReadVariableOpAdam/m/Adv/kernel*
_output_shapes

:2*
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
Ж
Adam/v/dense_9/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*&
shared_nameAdam/v/dense_9/kernel

)Adam/v/dense_9/kernel/Read/ReadVariableOpReadVariableOpAdam/v/dense_9/kernel*
_output_shapes

:22*
dtype0
Ж
Adam/m/dense_9/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*&
shared_nameAdam/m/dense_9/kernel

)Adam/m/dense_9/kernel/Read/ReadVariableOpReadVariableOpAdam/m/dense_9/kernel*
_output_shapes

:22*
dtype0
~
Adam/v/dense_8/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*$
shared_nameAdam/v/dense_8/bias
w
'Adam/v/dense_8/bias/Read/ReadVariableOpReadVariableOpAdam/v/dense_8/bias*
_output_shapes
:2*
dtype0
~
Adam/m/dense_8/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*$
shared_nameAdam/m/dense_8/bias
w
'Adam/m/dense_8/bias/Read/ReadVariableOpReadVariableOpAdam/m/dense_8/bias*
_output_shapes
:2*
dtype0
Ж
Adam/v/dense_8/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*&
shared_nameAdam/v/dense_8/kernel

)Adam/v/dense_8/kernel/Read/ReadVariableOpReadVariableOpAdam/v/dense_8/kernel*
_output_shapes

:22*
dtype0
Ж
Adam/m/dense_8/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*&
shared_nameAdam/m/dense_8/kernel

)Adam/m/dense_8/kernel/Read/ReadVariableOpReadVariableOpAdam/m/dense_8/kernel*
_output_shapes

:22*
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
Ж
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
Ж
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
Ж
Adam/v/dense_6/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*&
shared_nameAdam/v/dense_6/kernel

)Adam/v/dense_6/kernel/Read/ReadVariableOpReadVariableOpAdam/v/dense_6/kernel*
_output_shapes

:22*
dtype0
Ж
Adam/m/dense_6/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*&
shared_nameAdam/m/dense_6/kernel

)Adam/m/dense_6/kernel/Read/ReadVariableOpReadVariableOpAdam/m/dense_6/kernel*
_output_shapes

:22*
dtype0
~
Adam/v/dense_5/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*$
shared_nameAdam/v/dense_5/bias
w
'Adam/v/dense_5/bias/Read/ReadVariableOpReadVariableOpAdam/v/dense_5/bias*
_output_shapes
:2*
dtype0
~
Adam/m/dense_5/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*$
shared_nameAdam/m/dense_5/bias
w
'Adam/m/dense_5/bias/Read/ReadVariableOpReadVariableOpAdam/m/dense_5/bias*
_output_shapes
:2*
dtype0
Ж
Adam/v/dense_5/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*&
shared_nameAdam/v/dense_5/kernel

)Adam/v/dense_5/kernel/Read/ReadVariableOpReadVariableOpAdam/v/dense_5/kernel*
_output_shapes

:2*
dtype0
Ж
Adam/m/dense_5/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*&
shared_nameAdam/m/dense_5/kernel

)Adam/m/dense_5/kernel/Read/ReadVariableOpReadVariableOpAdam/m/dense_5/kernel*
_output_shapes

:2*
dtype0
v
Adam/v/Clf/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_nameAdam/v/Clf/bias
o
#Adam/v/Clf/bias/Read/ReadVariableOpReadVariableOpAdam/v/Clf/bias*
_output_shapes
:*
dtype0
v
Adam/m/Clf/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_nameAdam/m/Clf/bias
o
#Adam/m/Clf/bias/Read/ReadVariableOpReadVariableOpAdam/m/Clf/bias*
_output_shapes
:*
dtype0
~
Adam/v/Clf/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*"
shared_nameAdam/v/Clf/kernel
w
%Adam/v/Clf/kernel/Read/ReadVariableOpReadVariableOpAdam/v/Clf/kernel*
_output_shapes

:2*
dtype0
~
Adam/m/Clf/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*"
shared_nameAdam/m/Clf/kernel
w
%Adam/m/Clf/kernel/Read/ReadVariableOpReadVariableOpAdam/m/Clf/kernel*
_output_shapes

:2*
dtype0
~
Adam/v/dense_4/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*$
shared_nameAdam/v/dense_4/bias
w
'Adam/v/dense_4/bias/Read/ReadVariableOpReadVariableOpAdam/v/dense_4/bias*
_output_shapes
:2*
dtype0
~
Adam/m/dense_4/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*$
shared_nameAdam/m/dense_4/bias
w
'Adam/m/dense_4/bias/Read/ReadVariableOpReadVariableOpAdam/m/dense_4/bias*
_output_shapes
:2*
dtype0
Ж
Adam/v/dense_4/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*&
shared_nameAdam/v/dense_4/kernel

)Adam/v/dense_4/kernel/Read/ReadVariableOpReadVariableOpAdam/v/dense_4/kernel*
_output_shapes

:22*
dtype0
Ж
Adam/m/dense_4/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*&
shared_nameAdam/m/dense_4/kernel

)Adam/m/dense_4/kernel/Read/ReadVariableOpReadVariableOpAdam/m/dense_4/kernel*
_output_shapes

:22*
dtype0
~
Adam/v/dense_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*$
shared_nameAdam/v/dense_3/bias
w
'Adam/v/dense_3/bias/Read/ReadVariableOpReadVariableOpAdam/v/dense_3/bias*
_output_shapes
:2*
dtype0
~
Adam/m/dense_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*$
shared_nameAdam/m/dense_3/bias
w
'Adam/m/dense_3/bias/Read/ReadVariableOpReadVariableOpAdam/m/dense_3/bias*
_output_shapes
:2*
dtype0
Ж
Adam/v/dense_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:	2*&
shared_nameAdam/v/dense_3/kernel

)Adam/v/dense_3/kernel/Read/ReadVariableOpReadVariableOpAdam/v/dense_3/kernel*
_output_shapes

:	2*
dtype0
Ж
Adam/m/dense_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:	2*&
shared_nameAdam/m/dense_3/kernel

)Adam/m/dense_3/kernel/Read/ReadVariableOpReadVariableOpAdam/m/dense_3/kernel*
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
h
Adv/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_name
Adv/bias
a
Adv/bias/Read/ReadVariableOpReadVariableOpAdv/bias*
_output_shapes
:*
dtype0
p

Adv/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*
shared_name
Adv/kernel
i
Adv/kernel/Read/ReadVariableOpReadVariableOp
Adv/kernel*
_output_shapes

:2*
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
:22*
shared_namedense_9/kernel
q
"dense_9/kernel/Read/ReadVariableOpReadVariableOpdense_9/kernel*
_output_shapes

:22*
dtype0
p
dense_8/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*
shared_namedense_8/bias
i
 dense_8/bias/Read/ReadVariableOpReadVariableOpdense_8/bias*
_output_shapes
:2*
dtype0
x
dense_8/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*
shared_namedense_8/kernel
q
"dense_8/kernel/Read/ReadVariableOpReadVariableOpdense_8/kernel*
_output_shapes

:22*
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
:22*
shared_namedense_6/kernel
q
"dense_6/kernel/Read/ReadVariableOpReadVariableOpdense_6/kernel*
_output_shapes

:22*
dtype0
p
dense_5/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*
shared_namedense_5/bias
i
 dense_5/bias/Read/ReadVariableOpReadVariableOpdense_5/bias*
_output_shapes
:2*
dtype0
x
dense_5/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*
shared_namedense_5/kernel
q
"dense_5/kernel/Read/ReadVariableOpReadVariableOpdense_5/kernel*
_output_shapes

:2*
dtype0
h
Clf/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_name
Clf/bias
a
Clf/bias/Read/ReadVariableOpReadVariableOpClf/bias*
_output_shapes
:*
dtype0
p

Clf/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*
shared_name
Clf/kernel
i
Clf/kernel/Read/ReadVariableOpReadVariableOp
Clf/kernel*
_output_shapes

:2*
dtype0
p
dense_4/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*
shared_namedense_4/bias
i
 dense_4/bias/Read/ReadVariableOpReadVariableOpdense_4/bias*
_output_shapes
:2*
dtype0
x
dense_4/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:22*
shared_namedense_4/kernel
q
"dense_4/kernel/Read/ReadVariableOpReadVariableOpdense_4/kernel*
_output_shapes

:22*
dtype0
p
dense_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*
shared_namedense_3/bias
i
 dense_3/bias/Read/ReadVariableOpReadVariableOpdense_3/bias*
_output_shapes
:2*
dtype0
x
dense_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:	2*
shared_namedense_3/kernel
q
"dense_3/kernel/Read/ReadVariableOpReadVariableOpdense_3/kernel*
_output_shapes

:	2*
dtype0
z
serving_default_input_1Placeholder*'
_output_shapes
:€€€€€€€€€	*
dtype0*
shape:€€€€€€€€€	
z
serving_default_input_2Placeholder*'
_output_shapes
:€€€€€€€€€*
dtype0*
shape:€€€€€€€€€
И
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_1serving_default_input_2dense_3/kerneldense_3/biasdense_4/kerneldense_4/bias
Clf/kernelClf/biasdense_5/kerneldense_5/biasdense_6/kerneldense_6/biasdense_7/kerneldense_7/biasdense_8/kerneldense_8/biasdense_9/kerneldense_9/bias
Adv/kernelAdv/bias*
Tin
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:€€€€€€€€€:€€€€€€€€€*4
_read_only_resource_inputs
	
*2
config_proto" 

CPU

GPU2 *0J 8В *.
f)R'
%__inference_signature_wrapper_5670972

NoOpNoOp
Бu
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*Љt
value≤tBѓt B®t
м
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer-4
layer-5
layer-6
layer_with_weights-3
layer-7
	layer_with_weights-4
	layer-8

layer_with_weights-5

layer-9
layer_with_weights-6
layer-10
layer_with_weights-7
layer-11
layer_with_weights-8
layer-12
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses
_default_save_signature
	optimizer
loss

signatures*
* 
¶
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses

kernel
bias*
¶
 	variables
!trainable_variables
"regularization_losses
#	keras_api
$__call__
*%&call_and_return_all_conditional_losses

&kernel
'bias*
¶
(	variables
)trainable_variables
*regularization_losses
+	keras_api
,__call__
*-&call_and_return_all_conditional_losses

.kernel
/bias*
О
0	variables
1trainable_variables
2regularization_losses
3	keras_api
4__call__
*5&call_and_return_all_conditional_losses* 
* 
О
6	variables
7trainable_variables
8regularization_losses
9	keras_api
:__call__
*;&call_and_return_all_conditional_losses* 
¶
<	variables
=trainable_variables
>regularization_losses
?	keras_api
@__call__
*A&call_and_return_all_conditional_losses

Bkernel
Cbias*
¶
D	variables
Etrainable_variables
Fregularization_losses
G	keras_api
H__call__
*I&call_and_return_all_conditional_losses

Jkernel
Kbias*
¶
L	variables
Mtrainable_variables
Nregularization_losses
O	keras_api
P__call__
*Q&call_and_return_all_conditional_losses

Rkernel
Sbias*
¶
T	variables
Utrainable_variables
Vregularization_losses
W	keras_api
X__call__
*Y&call_and_return_all_conditional_losses

Zkernel
[bias*
¶
\	variables
]trainable_variables
^regularization_losses
_	keras_api
`__call__
*a&call_and_return_all_conditional_losses

bkernel
cbias*
¶
d	variables
etrainable_variables
fregularization_losses
g	keras_api
h__call__
*i&call_and_return_all_conditional_losses

jkernel
kbias*
К
0
1
&2
'3
.4
/5
B6
C7
J8
K9
R10
S11
Z12
[13
b14
c15
j16
k17*
К
0
1
&2
'3
.4
/5
B6
C7
J8
K9
R10
S11
Z12
[13
b14
c15
j16
k17*
* 
∞
lnon_trainable_variables

mlayers
nmetrics
olayer_regularization_losses
player_metrics
	variables
trainable_variables
regularization_losses
__call__
_default_save_signature
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses*
6
qtrace_0
rtrace_1
strace_2
ttrace_3* 
6
utrace_0
vtrace_1
wtrace_2
xtrace_3* 
* 
Б
y
_variables
z_iterations
{_learning_rate
|_index_dict
}
_momentums
~_velocities
_update_step_xla*
* 

Аserving_default* 

0
1*

0
1*
* 
Ш
Бnon_trainable_variables
Вlayers
Гmetrics
 Дlayer_regularization_losses
Еlayer_metrics
	variables
trainable_variables
regularization_losses
__call__
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses*

Жtrace_0* 

Зtrace_0* 
^X
VARIABLE_VALUEdense_3/kernel6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUEdense_3/bias4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE*

&0
'1*

&0
'1*
* 
Ш
Иnon_trainable_variables
Йlayers
Кmetrics
 Лlayer_regularization_losses
Мlayer_metrics
 	variables
!trainable_variables
"regularization_losses
$__call__
*%&call_and_return_all_conditional_losses
&%"call_and_return_conditional_losses*

Нtrace_0* 

Оtrace_0* 
^X
VARIABLE_VALUEdense_4/kernel6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUEdense_4/bias4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE*

.0
/1*

.0
/1*
* 
Ш
Пnon_trainable_variables
Рlayers
Сmetrics
 Тlayer_regularization_losses
Уlayer_metrics
(	variables
)trainable_variables
*regularization_losses
,__call__
*-&call_and_return_all_conditional_losses
&-"call_and_return_conditional_losses*

Фtrace_0* 

Хtrace_0* 
ZT
VARIABLE_VALUE
Clf/kernel6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE*
VP
VARIABLE_VALUEClf/bias4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUE*
* 
* 
* 
Ц
Цnon_trainable_variables
Чlayers
Шmetrics
 Щlayer_regularization_losses
Ъlayer_metrics
0	variables
1trainable_variables
2regularization_losses
4__call__
*5&call_and_return_all_conditional_losses
&5"call_and_return_conditional_losses* 

Ыtrace_0* 

Ьtrace_0* 
* 
* 
* 
Ц
Эnon_trainable_variables
Юlayers
Яmetrics
 †layer_regularization_losses
°layer_metrics
6	variables
7trainable_variables
8regularization_losses
:__call__
*;&call_and_return_all_conditional_losses
&;"call_and_return_conditional_losses* 

Ґtrace_0* 

£trace_0* 

B0
C1*

B0
C1*
* 
Ш
§non_trainable_variables
•layers
¶metrics
 Іlayer_regularization_losses
®layer_metrics
<	variables
=trainable_variables
>regularization_losses
@__call__
*A&call_and_return_all_conditional_losses
&A"call_and_return_conditional_losses*

©trace_0* 

™trace_0* 
^X
VARIABLE_VALUEdense_5/kernel6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUEdense_5/bias4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUE*

J0
K1*

J0
K1*
* 
Ш
Ђnon_trainable_variables
ђlayers
≠metrics
 Ѓlayer_regularization_losses
ѓlayer_metrics
D	variables
Etrainable_variables
Fregularization_losses
H__call__
*I&call_and_return_all_conditional_losses
&I"call_and_return_conditional_losses*

∞trace_0* 

±trace_0* 
^X
VARIABLE_VALUEdense_6/kernel6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUEdense_6/bias4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUE*

R0
S1*

R0
S1*
* 
Ш
≤non_trainable_variables
≥layers
іmetrics
 µlayer_regularization_losses
ґlayer_metrics
L	variables
Mtrainable_variables
Nregularization_losses
P__call__
*Q&call_and_return_all_conditional_losses
&Q"call_and_return_conditional_losses*

Јtrace_0* 

Єtrace_0* 
^X
VARIABLE_VALUEdense_7/kernel6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUEdense_7/bias4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUE*

Z0
[1*

Z0
[1*
* 
Ш
єnon_trainable_variables
Їlayers
їmetrics
 Љlayer_regularization_losses
љlayer_metrics
T	variables
Utrainable_variables
Vregularization_losses
X__call__
*Y&call_and_return_all_conditional_losses
&Y"call_and_return_conditional_losses*

Њtrace_0* 

њtrace_0* 
^X
VARIABLE_VALUEdense_8/kernel6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUEdense_8/bias4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUE*

b0
c1*

b0
c1*
* 
Ш
јnon_trainable_variables
Ѕlayers
¬metrics
 √layer_regularization_losses
ƒlayer_metrics
\	variables
]trainable_variables
^regularization_losses
`__call__
*a&call_and_return_all_conditional_losses
&a"call_and_return_conditional_losses*

≈trace_0* 

∆trace_0* 
^X
VARIABLE_VALUEdense_9/kernel6layer_with_weights-7/kernel/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUEdense_9/bias4layer_with_weights-7/bias/.ATTRIBUTES/VARIABLE_VALUE*

j0
k1*

j0
k1*
* 
Ш
«non_trainable_variables
»layers
…metrics
  layer_regularization_losses
Ћlayer_metrics
d	variables
etrainable_variables
fregularization_losses
h__call__
*i&call_and_return_all_conditional_losses
&i"call_and_return_conditional_losses*

ћtrace_0* 

Ќtrace_0* 
ZT
VARIABLE_VALUE
Adv/kernel6layer_with_weights-8/kernel/.ATTRIBUTES/VARIABLE_VALUE*
VP
VARIABLE_VALUEAdv/bias4layer_with_weights-8/bias/.ATTRIBUTES/VARIABLE_VALUE*
* 
b
0
1
2
3
4
5
6
7
	8

9
10
11
12*

ќ0
ѕ1
–2*
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
∆
z0
—1
“2
”3
‘4
’5
÷6
„7
Ў8
ў9
Џ10
џ11
№12
Ё13
ё14
я15
а16
б17
в18
г19
д20
е21
ж22
з23
и24
й25
к26
л27
м28
н29
о30
п31
р32
с33
т34
у35
ф36*
SM
VARIABLE_VALUE	iteration0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUElearning_rate3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
* 
Ь
—0
”1
’2
„3
ў4
џ5
Ё6
я7
б8
г9
е10
з11
й12
л13
н14
п15
с16
у17*
Ь
“0
‘1
÷2
Ў3
Џ4
№5
ё6
а7
в8
д9
ж10
и11
к12
м13
о14
р15
т16
ф17*
Ж
хtrace_0
цtrace_1
чtrace_2
шtrace_3
щtrace_4
ъtrace_5
ыtrace_6
ьtrace_7
эtrace_8
юtrace_9
€trace_10
Аtrace_11
Бtrace_12
Вtrace_13
Гtrace_14
Дtrace_15
Еtrace_16
Жtrace_17* 
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
<
З	variables
И	keras_api

Йtotal

Кcount*
<
Л	variables
М	keras_api

Нtotal

Оcount*
<
П	variables
Р	keras_api

Сtotal

Тcount*
`Z
VARIABLE_VALUEAdam/m/dense_3/kernel1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUEAdam/v/dense_3/kernel1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEAdam/m/dense_3/bias1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEAdam/v/dense_3/bias1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUEAdam/m/dense_4/kernel1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUE*
`Z
VARIABLE_VALUEAdam/v/dense_4/kernel1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEAdam/m/dense_4/bias1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEAdam/v/dense_4/bias1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUE*
\V
VARIABLE_VALUEAdam/m/Clf/kernel1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUE*
]W
VARIABLE_VALUEAdam/v/Clf/kernel2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/m/Clf/bias2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/Clf/bias2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/m/dense_5/kernel2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/v/dense_5/kernel2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/m/dense_5/bias2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/v/dense_5/bias2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/m/dense_6/kernel2optimizer/_variables/17/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/v/dense_6/kernel2optimizer/_variables/18/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/m/dense_6/bias2optimizer/_variables/19/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/v/dense_6/bias2optimizer/_variables/20/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/m/dense_7/kernel2optimizer/_variables/21/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/v/dense_7/kernel2optimizer/_variables/22/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/m/dense_7/bias2optimizer/_variables/23/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/v/dense_7/bias2optimizer/_variables/24/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/m/dense_8/kernel2optimizer/_variables/25/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/v/dense_8/kernel2optimizer/_variables/26/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/m/dense_8/bias2optimizer/_variables/27/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/v/dense_8/bias2optimizer/_variables/28/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/m/dense_9/kernel2optimizer/_variables/29/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/v/dense_9/kernel2optimizer/_variables/30/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/m/dense_9/bias2optimizer/_variables/31/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/v/dense_9/bias2optimizer/_variables/32/.ATTRIBUTES/VARIABLE_VALUE*
]W
VARIABLE_VALUEAdam/m/Adv/kernel2optimizer/_variables/33/.ATTRIBUTES/VARIABLE_VALUE*
]W
VARIABLE_VALUEAdam/v/Adv/kernel2optimizer/_variables/34/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/m/Adv/bias2optimizer/_variables/35/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/Adv/bias2optimizer/_variables/36/.ATTRIBUTES/VARIABLE_VALUE*
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
Й0
К1*

З	variables*
UO
VARIABLE_VALUEtotal_24keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEcount_24keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE*

Н0
О1*

Л	variables*
UO
VARIABLE_VALUEtotal_14keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEcount_14keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUE*

С0
Т1*

П	variables*
SM
VARIABLE_VALUEtotal4keras_api/metrics/2/total/.ATTRIBUTES/VARIABLE_VALUE*
SM
VARIABLE_VALUEcount4keras_api/metrics/2/count/.ATTRIBUTES/VARIABLE_VALUE*
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
щ
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename"dense_3/kernel/Read/ReadVariableOp dense_3/bias/Read/ReadVariableOp"dense_4/kernel/Read/ReadVariableOp dense_4/bias/Read/ReadVariableOpClf/kernel/Read/ReadVariableOpClf/bias/Read/ReadVariableOp"dense_5/kernel/Read/ReadVariableOp dense_5/bias/Read/ReadVariableOp"dense_6/kernel/Read/ReadVariableOp dense_6/bias/Read/ReadVariableOp"dense_7/kernel/Read/ReadVariableOp dense_7/bias/Read/ReadVariableOp"dense_8/kernel/Read/ReadVariableOp dense_8/bias/Read/ReadVariableOp"dense_9/kernel/Read/ReadVariableOp dense_9/bias/Read/ReadVariableOpAdv/kernel/Read/ReadVariableOpAdv/bias/Read/ReadVariableOpiteration/Read/ReadVariableOp!learning_rate/Read/ReadVariableOp)Adam/m/dense_3/kernel/Read/ReadVariableOp)Adam/v/dense_3/kernel/Read/ReadVariableOp'Adam/m/dense_3/bias/Read/ReadVariableOp'Adam/v/dense_3/bias/Read/ReadVariableOp)Adam/m/dense_4/kernel/Read/ReadVariableOp)Adam/v/dense_4/kernel/Read/ReadVariableOp'Adam/m/dense_4/bias/Read/ReadVariableOp'Adam/v/dense_4/bias/Read/ReadVariableOp%Adam/m/Clf/kernel/Read/ReadVariableOp%Adam/v/Clf/kernel/Read/ReadVariableOp#Adam/m/Clf/bias/Read/ReadVariableOp#Adam/v/Clf/bias/Read/ReadVariableOp)Adam/m/dense_5/kernel/Read/ReadVariableOp)Adam/v/dense_5/kernel/Read/ReadVariableOp'Adam/m/dense_5/bias/Read/ReadVariableOp'Adam/v/dense_5/bias/Read/ReadVariableOp)Adam/m/dense_6/kernel/Read/ReadVariableOp)Adam/v/dense_6/kernel/Read/ReadVariableOp'Adam/m/dense_6/bias/Read/ReadVariableOp'Adam/v/dense_6/bias/Read/ReadVariableOp)Adam/m/dense_7/kernel/Read/ReadVariableOp)Adam/v/dense_7/kernel/Read/ReadVariableOp'Adam/m/dense_7/bias/Read/ReadVariableOp'Adam/v/dense_7/bias/Read/ReadVariableOp)Adam/m/dense_8/kernel/Read/ReadVariableOp)Adam/v/dense_8/kernel/Read/ReadVariableOp'Adam/m/dense_8/bias/Read/ReadVariableOp'Adam/v/dense_8/bias/Read/ReadVariableOp)Adam/m/dense_9/kernel/Read/ReadVariableOp)Adam/v/dense_9/kernel/Read/ReadVariableOp'Adam/m/dense_9/bias/Read/ReadVariableOp'Adam/v/dense_9/bias/Read/ReadVariableOp%Adam/m/Adv/kernel/Read/ReadVariableOp%Adam/v/Adv/kernel/Read/ReadVariableOp#Adam/m/Adv/bias/Read/ReadVariableOp#Adam/v/Adv/bias/Read/ReadVariableOptotal_2/Read/ReadVariableOpcount_2/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOpConst*K
TinD
B2@	*
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
 __inference__traced_save_5671765
Ь
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenamedense_3/kerneldense_3/biasdense_4/kerneldense_4/bias
Clf/kernelClf/biasdense_5/kerneldense_5/biasdense_6/kerneldense_6/biasdense_7/kerneldense_7/biasdense_8/kerneldense_8/biasdense_9/kerneldense_9/bias
Adv/kernelAdv/bias	iterationlearning_rateAdam/m/dense_3/kernelAdam/v/dense_3/kernelAdam/m/dense_3/biasAdam/v/dense_3/biasAdam/m/dense_4/kernelAdam/v/dense_4/kernelAdam/m/dense_4/biasAdam/v/dense_4/biasAdam/m/Clf/kernelAdam/v/Clf/kernelAdam/m/Clf/biasAdam/v/Clf/biasAdam/m/dense_5/kernelAdam/v/dense_5/kernelAdam/m/dense_5/biasAdam/v/dense_5/biasAdam/m/dense_6/kernelAdam/v/dense_6/kernelAdam/m/dense_6/biasAdam/v/dense_6/biasAdam/m/dense_7/kernelAdam/v/dense_7/kernelAdam/m/dense_7/biasAdam/v/dense_7/biasAdam/m/dense_8/kernelAdam/v/dense_8/kernelAdam/m/dense_8/biasAdam/v/dense_8/biasAdam/m/dense_9/kernelAdam/v/dense_9/kernelAdam/m/dense_9/biasAdam/v/dense_9/biasAdam/m/Adv/kernelAdam/v/Adv/kernelAdam/m/Adv/biasAdam/v/Adv/biastotal_2count_2total_1count_1totalcount*J
TinC
A2?*
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
#__inference__traced_restore_5671961ёэ	
Ы

х
D__inference_dense_8_layer_call_and_return_conditional_losses_5670447

inputs0
matmul_readvariableop_resource:22-
biasadd_readvariableop_resource:2
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:22*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
і
н
)__inference_model_1_layer_call_fn_5671016
inputs_0
inputs_1
unknown:	2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
	unknown_5:2
	unknown_6:2
	unknown_7:22
	unknown_8:2
	unknown_9:22

unknown_10:2

unknown_11:22

unknown_12:2

unknown_13:22

unknown_14:2

unknown_15:2

unknown_16:
identity

identity_1ИҐStatefulPartitionedCall÷
StatefulPartitionedCallStatefulPartitionedCallinputs_0inputs_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:€€€€€€€€€:€€€€€€€€€*4
_read_only_resource_inputs
	
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_model_1_layer_call_and_return_conditional_losses_5670489o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€q

Identity_1Identity StatefulPartitionedCall:output:1^NoOp*
T0*'
_output_shapes
:€€€€€€€€€`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*]
_input_shapesL
J:€€€€€€€€€	:€€€€€€€€€: : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:Q M
'
_output_shapes
:€€€€€€€€€	
"
_user_specified_name
inputs_0:QM
'
_output_shapes
:€€€€€€€€€
"
_user_specified_name
inputs_1
ђ
\
$__inference_internal_grad_fn_5671665
result_grads_0
result_grads_1
identityL
NegNegresult_grads_0*
T0*'
_output_shapes
:€€€€€€€€€J
mul/yConst*
_output_shapes
: *
dtype0*
valueB
 *   AU
mulMulNeg:y:0mul/y:output:0*
T0*'
_output_shapes
:€€€€€€€€€O
IdentityIdentitymul:z:0*
T0*'
_output_shapes
:€€€€€€€€€"
identityIdentity:output:0*9
_input_shapes(
&:€€€€€€€€€:€€€€€€€€€:W S
'
_output_shapes
:€€€€€€€€€
(
_user_specified_nameresult_grads_0:WS
'
_output_shapes
:€€€€€€€€€
(
_user_specified_nameresult_grads_1
Ц

с
@__inference_Clf_layer_call_and_return_conditional_losses_5670359

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:2*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€V
SigmoidSigmoidBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€Z
IdentityIdentitySigmoid:y:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
ђ€
–$
#__inference__traced_restore_5671961
file_prefix1
assignvariableop_dense_3_kernel:	2-
assignvariableop_1_dense_3_bias:23
!assignvariableop_2_dense_4_kernel:22-
assignvariableop_3_dense_4_bias:2/
assignvariableop_4_clf_kernel:2)
assignvariableop_5_clf_bias:3
!assignvariableop_6_dense_5_kernel:2-
assignvariableop_7_dense_5_bias:23
!assignvariableop_8_dense_6_kernel:22-
assignvariableop_9_dense_6_bias:24
"assignvariableop_10_dense_7_kernel:22.
 assignvariableop_11_dense_7_bias:24
"assignvariableop_12_dense_8_kernel:22.
 assignvariableop_13_dense_8_bias:24
"assignvariableop_14_dense_9_kernel:22.
 assignvariableop_15_dense_9_bias:20
assignvariableop_16_adv_kernel:2*
assignvariableop_17_adv_bias:'
assignvariableop_18_iteration:	 +
!assignvariableop_19_learning_rate: ;
)assignvariableop_20_adam_m_dense_3_kernel:	2;
)assignvariableop_21_adam_v_dense_3_kernel:	25
'assignvariableop_22_adam_m_dense_3_bias:25
'assignvariableop_23_adam_v_dense_3_bias:2;
)assignvariableop_24_adam_m_dense_4_kernel:22;
)assignvariableop_25_adam_v_dense_4_kernel:225
'assignvariableop_26_adam_m_dense_4_bias:25
'assignvariableop_27_adam_v_dense_4_bias:27
%assignvariableop_28_adam_m_clf_kernel:27
%assignvariableop_29_adam_v_clf_kernel:21
#assignvariableop_30_adam_m_clf_bias:1
#assignvariableop_31_adam_v_clf_bias:;
)assignvariableop_32_adam_m_dense_5_kernel:2;
)assignvariableop_33_adam_v_dense_5_kernel:25
'assignvariableop_34_adam_m_dense_5_bias:25
'assignvariableop_35_adam_v_dense_5_bias:2;
)assignvariableop_36_adam_m_dense_6_kernel:22;
)assignvariableop_37_adam_v_dense_6_kernel:225
'assignvariableop_38_adam_m_dense_6_bias:25
'assignvariableop_39_adam_v_dense_6_bias:2;
)assignvariableop_40_adam_m_dense_7_kernel:22;
)assignvariableop_41_adam_v_dense_7_kernel:225
'assignvariableop_42_adam_m_dense_7_bias:25
'assignvariableop_43_adam_v_dense_7_bias:2;
)assignvariableop_44_adam_m_dense_8_kernel:22;
)assignvariableop_45_adam_v_dense_8_kernel:225
'assignvariableop_46_adam_m_dense_8_bias:25
'assignvariableop_47_adam_v_dense_8_bias:2;
)assignvariableop_48_adam_m_dense_9_kernel:22;
)assignvariableop_49_adam_v_dense_9_kernel:225
'assignvariableop_50_adam_m_dense_9_bias:25
'assignvariableop_51_adam_v_dense_9_bias:27
%assignvariableop_52_adam_m_adv_kernel:27
%assignvariableop_53_adam_v_adv_kernel:21
#assignvariableop_54_adam_m_adv_bias:1
#assignvariableop_55_adam_v_adv_bias:%
assignvariableop_56_total_2: %
assignvariableop_57_count_2: %
assignvariableop_58_total_1: %
assignvariableop_59_count_1: #
assignvariableop_60_total: #
assignvariableop_61_count: 
identity_63ИҐAssignVariableOpҐAssignVariableOp_1ҐAssignVariableOp_10ҐAssignVariableOp_11ҐAssignVariableOp_12ҐAssignVariableOp_13ҐAssignVariableOp_14ҐAssignVariableOp_15ҐAssignVariableOp_16ҐAssignVariableOp_17ҐAssignVariableOp_18ҐAssignVariableOp_19ҐAssignVariableOp_2ҐAssignVariableOp_20ҐAssignVariableOp_21ҐAssignVariableOp_22ҐAssignVariableOp_23ҐAssignVariableOp_24ҐAssignVariableOp_25ҐAssignVariableOp_26ҐAssignVariableOp_27ҐAssignVariableOp_28ҐAssignVariableOp_29ҐAssignVariableOp_3ҐAssignVariableOp_30ҐAssignVariableOp_31ҐAssignVariableOp_32ҐAssignVariableOp_33ҐAssignVariableOp_34ҐAssignVariableOp_35ҐAssignVariableOp_36ҐAssignVariableOp_37ҐAssignVariableOp_38ҐAssignVariableOp_39ҐAssignVariableOp_4ҐAssignVariableOp_40ҐAssignVariableOp_41ҐAssignVariableOp_42ҐAssignVariableOp_43ҐAssignVariableOp_44ҐAssignVariableOp_45ҐAssignVariableOp_46ҐAssignVariableOp_47ҐAssignVariableOp_48ҐAssignVariableOp_49ҐAssignVariableOp_5ҐAssignVariableOp_50ҐAssignVariableOp_51ҐAssignVariableOp_52ҐAssignVariableOp_53ҐAssignVariableOp_54ҐAssignVariableOp_55ҐAssignVariableOp_56ҐAssignVariableOp_57ҐAssignVariableOp_58ҐAssignVariableOp_59ҐAssignVariableOp_6ҐAssignVariableOp_60ҐAssignVariableOp_61ҐAssignVariableOp_7ҐAssignVariableOp_8ҐAssignVariableOp_9Ё
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:?*
dtype0*Г
valueщBц?B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-7/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-7/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-8/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-8/bias/.ATTRIBUTES/VARIABLE_VALUEB0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUEB3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/17/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/18/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/19/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/20/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/21/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/22/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/23/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/24/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/25/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/26/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/27/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/28/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/29/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/30/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/31/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/32/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/33/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/34/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/35/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/36/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/2/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/2/count/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPHс
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:?*
dtype0*У
valueЙBЖ?B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B №
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*Т
_output_shapes€
ь:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*M
dtypesC
A2?	[
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:≤
AssignVariableOpAssignVariableOpassignvariableop_dense_3_kernelIdentity:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:ґ
AssignVariableOp_1AssignVariableOpassignvariableop_1_dense_3_biasIdentity_1:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:Є
AssignVariableOp_2AssignVariableOp!assignvariableop_2_dense_4_kernelIdentity_2:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:ґ
AssignVariableOp_3AssignVariableOpassignvariableop_3_dense_4_biasIdentity_3:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:і
AssignVariableOp_4AssignVariableOpassignvariableop_4_clf_kernelIdentity_4:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:≤
AssignVariableOp_5AssignVariableOpassignvariableop_5_clf_biasIdentity_5:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:Є
AssignVariableOp_6AssignVariableOp!assignvariableop_6_dense_5_kernelIdentity_6:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:ґ
AssignVariableOp_7AssignVariableOpassignvariableop_7_dense_5_biasIdentity_7:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:Є
AssignVariableOp_8AssignVariableOp!assignvariableop_8_dense_6_kernelIdentity_8:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:ґ
AssignVariableOp_9AssignVariableOpassignvariableop_9_dense_6_biasIdentity_9:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:ї
AssignVariableOp_10AssignVariableOp"assignvariableop_10_dense_7_kernelIdentity_10:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:є
AssignVariableOp_11AssignVariableOp assignvariableop_11_dense_7_biasIdentity_11:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:ї
AssignVariableOp_12AssignVariableOp"assignvariableop_12_dense_8_kernelIdentity_12:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:є
AssignVariableOp_13AssignVariableOp assignvariableop_13_dense_8_biasIdentity_13:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:ї
AssignVariableOp_14AssignVariableOp"assignvariableop_14_dense_9_kernelIdentity_14:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:є
AssignVariableOp_15AssignVariableOp assignvariableop_15_dense_9_biasIdentity_15:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:Ј
AssignVariableOp_16AssignVariableOpassignvariableop_16_adv_kernelIdentity_16:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:µ
AssignVariableOp_17AssignVariableOpassignvariableop_17_adv_biasIdentity_17:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0	*
_output_shapes
:ґ
AssignVariableOp_18AssignVariableOpassignvariableop_18_iterationIdentity_18:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0	_
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:Ї
AssignVariableOp_19AssignVariableOp!assignvariableop_19_learning_rateIdentity_19:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:¬
AssignVariableOp_20AssignVariableOp)assignvariableop_20_adam_m_dense_3_kernelIdentity_20:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:¬
AssignVariableOp_21AssignVariableOp)assignvariableop_21_adam_v_dense_3_kernelIdentity_21:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:ј
AssignVariableOp_22AssignVariableOp'assignvariableop_22_adam_m_dense_3_biasIdentity_22:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:ј
AssignVariableOp_23AssignVariableOp'assignvariableop_23_adam_v_dense_3_biasIdentity_23:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:¬
AssignVariableOp_24AssignVariableOp)assignvariableop_24_adam_m_dense_4_kernelIdentity_24:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:¬
AssignVariableOp_25AssignVariableOp)assignvariableop_25_adam_v_dense_4_kernelIdentity_25:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:ј
AssignVariableOp_26AssignVariableOp'assignvariableop_26_adam_m_dense_4_biasIdentity_26:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:ј
AssignVariableOp_27AssignVariableOp'assignvariableop_27_adam_v_dense_4_biasIdentity_27:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:Њ
AssignVariableOp_28AssignVariableOp%assignvariableop_28_adam_m_clf_kernelIdentity_28:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:Њ
AssignVariableOp_29AssignVariableOp%assignvariableop_29_adam_v_clf_kernelIdentity_29:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:Љ
AssignVariableOp_30AssignVariableOp#assignvariableop_30_adam_m_clf_biasIdentity_30:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:Љ
AssignVariableOp_31AssignVariableOp#assignvariableop_31_adam_v_clf_biasIdentity_31:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:¬
AssignVariableOp_32AssignVariableOp)assignvariableop_32_adam_m_dense_5_kernelIdentity_32:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:¬
AssignVariableOp_33AssignVariableOp)assignvariableop_33_adam_v_dense_5_kernelIdentity_33:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:ј
AssignVariableOp_34AssignVariableOp'assignvariableop_34_adam_m_dense_5_biasIdentity_34:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:ј
AssignVariableOp_35AssignVariableOp'assignvariableop_35_adam_v_dense_5_biasIdentity_35:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_36IdentityRestoreV2:tensors:36"/device:CPU:0*
T0*
_output_shapes
:¬
AssignVariableOp_36AssignVariableOp)assignvariableop_36_adam_m_dense_6_kernelIdentity_36:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_37IdentityRestoreV2:tensors:37"/device:CPU:0*
T0*
_output_shapes
:¬
AssignVariableOp_37AssignVariableOp)assignvariableop_37_adam_v_dense_6_kernelIdentity_37:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_38IdentityRestoreV2:tensors:38"/device:CPU:0*
T0*
_output_shapes
:ј
AssignVariableOp_38AssignVariableOp'assignvariableop_38_adam_m_dense_6_biasIdentity_38:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_39IdentityRestoreV2:tensors:39"/device:CPU:0*
T0*
_output_shapes
:ј
AssignVariableOp_39AssignVariableOp'assignvariableop_39_adam_v_dense_6_biasIdentity_39:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_40IdentityRestoreV2:tensors:40"/device:CPU:0*
T0*
_output_shapes
:¬
AssignVariableOp_40AssignVariableOp)assignvariableop_40_adam_m_dense_7_kernelIdentity_40:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_41IdentityRestoreV2:tensors:41"/device:CPU:0*
T0*
_output_shapes
:¬
AssignVariableOp_41AssignVariableOp)assignvariableop_41_adam_v_dense_7_kernelIdentity_41:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_42IdentityRestoreV2:tensors:42"/device:CPU:0*
T0*
_output_shapes
:ј
AssignVariableOp_42AssignVariableOp'assignvariableop_42_adam_m_dense_7_biasIdentity_42:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_43IdentityRestoreV2:tensors:43"/device:CPU:0*
T0*
_output_shapes
:ј
AssignVariableOp_43AssignVariableOp'assignvariableop_43_adam_v_dense_7_biasIdentity_43:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_44IdentityRestoreV2:tensors:44"/device:CPU:0*
T0*
_output_shapes
:¬
AssignVariableOp_44AssignVariableOp)assignvariableop_44_adam_m_dense_8_kernelIdentity_44:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_45IdentityRestoreV2:tensors:45"/device:CPU:0*
T0*
_output_shapes
:¬
AssignVariableOp_45AssignVariableOp)assignvariableop_45_adam_v_dense_8_kernelIdentity_45:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_46IdentityRestoreV2:tensors:46"/device:CPU:0*
T0*
_output_shapes
:ј
AssignVariableOp_46AssignVariableOp'assignvariableop_46_adam_m_dense_8_biasIdentity_46:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_47IdentityRestoreV2:tensors:47"/device:CPU:0*
T0*
_output_shapes
:ј
AssignVariableOp_47AssignVariableOp'assignvariableop_47_adam_v_dense_8_biasIdentity_47:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_48IdentityRestoreV2:tensors:48"/device:CPU:0*
T0*
_output_shapes
:¬
AssignVariableOp_48AssignVariableOp)assignvariableop_48_adam_m_dense_9_kernelIdentity_48:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_49IdentityRestoreV2:tensors:49"/device:CPU:0*
T0*
_output_shapes
:¬
AssignVariableOp_49AssignVariableOp)assignvariableop_49_adam_v_dense_9_kernelIdentity_49:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_50IdentityRestoreV2:tensors:50"/device:CPU:0*
T0*
_output_shapes
:ј
AssignVariableOp_50AssignVariableOp'assignvariableop_50_adam_m_dense_9_biasIdentity_50:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_51IdentityRestoreV2:tensors:51"/device:CPU:0*
T0*
_output_shapes
:ј
AssignVariableOp_51AssignVariableOp'assignvariableop_51_adam_v_dense_9_biasIdentity_51:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_52IdentityRestoreV2:tensors:52"/device:CPU:0*
T0*
_output_shapes
:Њ
AssignVariableOp_52AssignVariableOp%assignvariableop_52_adam_m_adv_kernelIdentity_52:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_53IdentityRestoreV2:tensors:53"/device:CPU:0*
T0*
_output_shapes
:Њ
AssignVariableOp_53AssignVariableOp%assignvariableop_53_adam_v_adv_kernelIdentity_53:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_54IdentityRestoreV2:tensors:54"/device:CPU:0*
T0*
_output_shapes
:Љ
AssignVariableOp_54AssignVariableOp#assignvariableop_54_adam_m_adv_biasIdentity_54:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_55IdentityRestoreV2:tensors:55"/device:CPU:0*
T0*
_output_shapes
:Љ
AssignVariableOp_55AssignVariableOp#assignvariableop_55_adam_v_adv_biasIdentity_55:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_56IdentityRestoreV2:tensors:56"/device:CPU:0*
T0*
_output_shapes
:і
AssignVariableOp_56AssignVariableOpassignvariableop_56_total_2Identity_56:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_57IdentityRestoreV2:tensors:57"/device:CPU:0*
T0*
_output_shapes
:і
AssignVariableOp_57AssignVariableOpassignvariableop_57_count_2Identity_57:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_58IdentityRestoreV2:tensors:58"/device:CPU:0*
T0*
_output_shapes
:і
AssignVariableOp_58AssignVariableOpassignvariableop_58_total_1Identity_58:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_59IdentityRestoreV2:tensors:59"/device:CPU:0*
T0*
_output_shapes
:і
AssignVariableOp_59AssignVariableOpassignvariableop_59_count_1Identity_59:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_60IdentityRestoreV2:tensors:60"/device:CPU:0*
T0*
_output_shapes
:≤
AssignVariableOp_60AssignVariableOpassignvariableop_60_totalIdentity_60:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_61IdentityRestoreV2:tensors:61"/device:CPU:0*
T0*
_output_shapes
:≤
AssignVariableOp_61AssignVariableOpassignvariableop_61_countIdentity_61:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0Y
NoOpNoOp"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 £
Identity_62Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_54^AssignVariableOp_55^AssignVariableOp_56^AssignVariableOp_57^AssignVariableOp_58^AssignVariableOp_59^AssignVariableOp_6^AssignVariableOp_60^AssignVariableOp_61^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: W
Identity_63IdentityIdentity_62:output:0^NoOp_1*
T0*
_output_shapes
: Р
NoOp_1NoOp^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_54^AssignVariableOp_55^AssignVariableOp_56^AssignVariableOp_57^AssignVariableOp_58^AssignVariableOp_59^AssignVariableOp_6^AssignVariableOp_60^AssignVariableOp_61^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9*"
_acd_function_control_output(*
_output_shapes
 "#
identity_63Identity_63:output:0*Т
_input_shapesА
~: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 2$
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
AssignVariableOp_23AssignVariableOp_232*
AssignVariableOp_24AssignVariableOp_242*
AssignVariableOp_25AssignVariableOp_252*
AssignVariableOp_26AssignVariableOp_262*
AssignVariableOp_27AssignVariableOp_272*
AssignVariableOp_28AssignVariableOp_282*
AssignVariableOp_29AssignVariableOp_292(
AssignVariableOp_3AssignVariableOp_32*
AssignVariableOp_30AssignVariableOp_302*
AssignVariableOp_31AssignVariableOp_312*
AssignVariableOp_32AssignVariableOp_322*
AssignVariableOp_33AssignVariableOp_332*
AssignVariableOp_34AssignVariableOp_342*
AssignVariableOp_35AssignVariableOp_352*
AssignVariableOp_36AssignVariableOp_362*
AssignVariableOp_37AssignVariableOp_372*
AssignVariableOp_38AssignVariableOp_382*
AssignVariableOp_39AssignVariableOp_392(
AssignVariableOp_4AssignVariableOp_42*
AssignVariableOp_40AssignVariableOp_402*
AssignVariableOp_41AssignVariableOp_412*
AssignVariableOp_42AssignVariableOp_422*
AssignVariableOp_43AssignVariableOp_432*
AssignVariableOp_44AssignVariableOp_442*
AssignVariableOp_45AssignVariableOp_452*
AssignVariableOp_46AssignVariableOp_462*
AssignVariableOp_47AssignVariableOp_472*
AssignVariableOp_48AssignVariableOp_482*
AssignVariableOp_49AssignVariableOp_492(
AssignVariableOp_5AssignVariableOp_52*
AssignVariableOp_50AssignVariableOp_502*
AssignVariableOp_51AssignVariableOp_512*
AssignVariableOp_52AssignVariableOp_522*
AssignVariableOp_53AssignVariableOp_532*
AssignVariableOp_54AssignVariableOp_542*
AssignVariableOp_55AssignVariableOp_552*
AssignVariableOp_56AssignVariableOp_562*
AssignVariableOp_57AssignVariableOp_572*
AssignVariableOp_58AssignVariableOp_582*
AssignVariableOp_59AssignVariableOp_592(
AssignVariableOp_6AssignVariableOp_62*
AssignVariableOp_60AssignVariableOp_602*
AssignVariableOp_61AssignVariableOp_612(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_9:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
є
P
$__inference__update_step_xla_5671257
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
њ
Т
%__inference_Adv_layer_call_fn_5671498

inputs
unknown:2
	unknown_0:
identityИҐStatefulPartitionedCallЏ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *I
fDRB
@__inference_Adv_layer_call_and_return_conditional_losses_5670481o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
є
P
$__inference__update_step_xla_5671227
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
Ц

с
@__inference_Clf_layer_call_and_return_conditional_losses_5671362

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:2*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€V
SigmoidSigmoidBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€Z
IdentityIdentitySigmoid:y:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
Ы

х
D__inference_dense_6_layer_call_and_return_conditional_losses_5671429

inputs0
matmul_readvariableop_resource:22-
biasadd_readvariableop_resource:2
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:22*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
є
P
$__inference__update_step_xla_5671277
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
Ы

х
D__inference_dense_8_layer_call_and_return_conditional_losses_5671469

inputs0
matmul_readvariableop_resource:22-
biasadd_readvariableop_resource:2
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:22*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
ђ
\
$__inference_internal_grad_fn_5671647
result_grads_0
result_grads_1
identityL
NegNegresult_grads_0*
T0*'
_output_shapes
:€€€€€€€€€J
mul/yConst*
_output_shapes
: *
dtype0*
valueB
 *   AU
mulMulNeg:y:0mul/y:output:0*
T0*'
_output_shapes
:€€€€€€€€€O
IdentityIdentitymul:z:0*
T0*'
_output_shapes
:€€€€€€€€€"
identityIdentity:output:0*9
_input_shapes(
&:€€€€€€€€€:€€€€€€€€€:W S
'
_output_shapes
:€€€€€€€€€
(
_user_specified_nameresult_grads_0:WS
'
_output_shapes
:€€€€€€€€€
(
_user_specified_nameresult_grads_1
Ы

х
D__inference_dense_9_layer_call_and_return_conditional_losses_5670464

inputs0
matmul_readvariableop_resource:22-
biasadd_readvariableop_resource:2
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:22*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
ћn
±
 __inference__traced_save_5671765
file_prefix-
)savev2_dense_3_kernel_read_readvariableop+
'savev2_dense_3_bias_read_readvariableop-
)savev2_dense_4_kernel_read_readvariableop+
'savev2_dense_4_bias_read_readvariableop)
%savev2_clf_kernel_read_readvariableop'
#savev2_clf_bias_read_readvariableop-
)savev2_dense_5_kernel_read_readvariableop+
'savev2_dense_5_bias_read_readvariableop-
)savev2_dense_6_kernel_read_readvariableop+
'savev2_dense_6_bias_read_readvariableop-
)savev2_dense_7_kernel_read_readvariableop+
'savev2_dense_7_bias_read_readvariableop-
)savev2_dense_8_kernel_read_readvariableop+
'savev2_dense_8_bias_read_readvariableop-
)savev2_dense_9_kernel_read_readvariableop+
'savev2_dense_9_bias_read_readvariableop)
%savev2_adv_kernel_read_readvariableop'
#savev2_adv_bias_read_readvariableop(
$savev2_iteration_read_readvariableop	,
(savev2_learning_rate_read_readvariableop4
0savev2_adam_m_dense_3_kernel_read_readvariableop4
0savev2_adam_v_dense_3_kernel_read_readvariableop2
.savev2_adam_m_dense_3_bias_read_readvariableop2
.savev2_adam_v_dense_3_bias_read_readvariableop4
0savev2_adam_m_dense_4_kernel_read_readvariableop4
0savev2_adam_v_dense_4_kernel_read_readvariableop2
.savev2_adam_m_dense_4_bias_read_readvariableop2
.savev2_adam_v_dense_4_bias_read_readvariableop0
,savev2_adam_m_clf_kernel_read_readvariableop0
,savev2_adam_v_clf_kernel_read_readvariableop.
*savev2_adam_m_clf_bias_read_readvariableop.
*savev2_adam_v_clf_bias_read_readvariableop4
0savev2_adam_m_dense_5_kernel_read_readvariableop4
0savev2_adam_v_dense_5_kernel_read_readvariableop2
.savev2_adam_m_dense_5_bias_read_readvariableop2
.savev2_adam_v_dense_5_bias_read_readvariableop4
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
.savev2_adam_v_dense_8_bias_read_readvariableop4
0savev2_adam_m_dense_9_kernel_read_readvariableop4
0savev2_adam_v_dense_9_kernel_read_readvariableop2
.savev2_adam_m_dense_9_bias_read_readvariableop2
.savev2_adam_v_dense_9_bias_read_readvariableop0
,savev2_adam_m_adv_kernel_read_readvariableop0
,savev2_adam_v_adv_kernel_read_readvariableop.
*savev2_adam_m_adv_bias_read_readvariableop.
*savev2_adam_v_adv_bias_read_readvariableop&
"savev2_total_2_read_readvariableop&
"savev2_count_2_read_readvariableop&
"savev2_total_1_read_readvariableop&
"savev2_count_1_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop
savev2_const

identity_1ИҐMergeV2Checkpointsw
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
: Џ
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:?*
dtype0*Г
valueщBц?B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-7/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-7/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-8/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-8/bias/.ATTRIBUTES/VARIABLE_VALUEB0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUEB3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/17/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/18/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/19/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/20/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/21/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/22/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/23/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/24/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/25/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/26/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/27/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/28/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/29/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/30/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/31/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/32/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/33/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/34/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/35/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/36/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/2/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/2/count/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPHо
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:?*
dtype0*У
valueЙBЖ?B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B в
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0)savev2_dense_3_kernel_read_readvariableop'savev2_dense_3_bias_read_readvariableop)savev2_dense_4_kernel_read_readvariableop'savev2_dense_4_bias_read_readvariableop%savev2_clf_kernel_read_readvariableop#savev2_clf_bias_read_readvariableop)savev2_dense_5_kernel_read_readvariableop'savev2_dense_5_bias_read_readvariableop)savev2_dense_6_kernel_read_readvariableop'savev2_dense_6_bias_read_readvariableop)savev2_dense_7_kernel_read_readvariableop'savev2_dense_7_bias_read_readvariableop)savev2_dense_8_kernel_read_readvariableop'savev2_dense_8_bias_read_readvariableop)savev2_dense_9_kernel_read_readvariableop'savev2_dense_9_bias_read_readvariableop%savev2_adv_kernel_read_readvariableop#savev2_adv_bias_read_readvariableop$savev2_iteration_read_readvariableop(savev2_learning_rate_read_readvariableop0savev2_adam_m_dense_3_kernel_read_readvariableop0savev2_adam_v_dense_3_kernel_read_readvariableop.savev2_adam_m_dense_3_bias_read_readvariableop.savev2_adam_v_dense_3_bias_read_readvariableop0savev2_adam_m_dense_4_kernel_read_readvariableop0savev2_adam_v_dense_4_kernel_read_readvariableop.savev2_adam_m_dense_4_bias_read_readvariableop.savev2_adam_v_dense_4_bias_read_readvariableop,savev2_adam_m_clf_kernel_read_readvariableop,savev2_adam_v_clf_kernel_read_readvariableop*savev2_adam_m_clf_bias_read_readvariableop*savev2_adam_v_clf_bias_read_readvariableop0savev2_adam_m_dense_5_kernel_read_readvariableop0savev2_adam_v_dense_5_kernel_read_readvariableop.savev2_adam_m_dense_5_bias_read_readvariableop.savev2_adam_v_dense_5_bias_read_readvariableop0savev2_adam_m_dense_6_kernel_read_readvariableop0savev2_adam_v_dense_6_kernel_read_readvariableop.savev2_adam_m_dense_6_bias_read_readvariableop.savev2_adam_v_dense_6_bias_read_readvariableop0savev2_adam_m_dense_7_kernel_read_readvariableop0savev2_adam_v_dense_7_kernel_read_readvariableop.savev2_adam_m_dense_7_bias_read_readvariableop.savev2_adam_v_dense_7_bias_read_readvariableop0savev2_adam_m_dense_8_kernel_read_readvariableop0savev2_adam_v_dense_8_kernel_read_readvariableop.savev2_adam_m_dense_8_bias_read_readvariableop.savev2_adam_v_dense_8_bias_read_readvariableop0savev2_adam_m_dense_9_kernel_read_readvariableop0savev2_adam_v_dense_9_kernel_read_readvariableop.savev2_adam_m_dense_9_bias_read_readvariableop.savev2_adam_v_dense_9_bias_read_readvariableop,savev2_adam_m_adv_kernel_read_readvariableop,savev2_adam_v_adv_kernel_read_readvariableop*savev2_adam_m_adv_bias_read_readvariableop*savev2_adam_v_adv_bias_read_readvariableop"savev2_total_2_read_readvariableop"savev2_count_2_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableopsavev2_const"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *M
dtypesC
A2?	Р
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:≥
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

identity_1Identity_1:output:0*ў
_input_shapes«
ƒ: :	2:2:22:2:2::2:2:22:2:22:2:22:2:22:2:2:: : :	2:	2:2:2:22:22:2:2:2:2:::2:2:2:2:22:22:2:2:22:22:2:2:22:22:2:2:22:22:2:2:2:2::: : : : : : : 2(
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
::$ 

_output_shapes

:2: 

_output_shapes
:2:$	 

_output_shapes

:22: 


_output_shapes
:2:$ 

_output_shapes

:22: 

_output_shapes
:2:$ 

_output_shapes

:22: 

_output_shapes
:2:$ 

_output_shapes

:22: 

_output_shapes
:2:$ 

_output_shapes

:2: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :$ 

_output_shapes

:	2:$ 

_output_shapes

:	2: 

_output_shapes
:2: 

_output_shapes
:2:$ 

_output_shapes

:22:$ 

_output_shapes

:22: 

_output_shapes
:2: 

_output_shapes
:2:$ 

_output_shapes

:2:$ 

_output_shapes

:2: 

_output_shapes
::  

_output_shapes
::$! 

_output_shapes

:2:$" 

_output_shapes

:2: #

_output_shapes
:2: $

_output_shapes
:2:$% 

_output_shapes

:22:$& 

_output_shapes

:22: '

_output_shapes
:2: (

_output_shapes
:2:$) 

_output_shapes

:22:$* 

_output_shapes

:22: +

_output_shapes
:2: ,

_output_shapes
:2:$- 

_output_shapes

:22:$. 

_output_shapes

:22: /

_output_shapes
:2: 0

_output_shapes
:2:$1 

_output_shapes

:22:$2 

_output_shapes

:22: 3

_output_shapes
:2: 4

_output_shapes
:2:$5 

_output_shapes

:2:$6 

_output_shapes

:2: 7

_output_shapes
:: 8

_output_shapes
::9

_output_shapes
: ::

_output_shapes
: :;

_output_shapes
: :<

_output_shapes
: :=

_output_shapes
: :>

_output_shapes
: :?

_output_shapes
: 
«
Ц
)__inference_dense_3_layer_call_fn_5671311

inputs
unknown:	2
	unknown_0:2
identityИҐStatefulPartitionedCallё
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_3_layer_call_and_return_conditional_losses_5670325o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€	: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:€€€€€€€€€	
 
_user_specified_nameinputs
“
b
I__inference_grad_reverse_layer_call_and_return_conditional_losses_5671376
x

identity_2I
IdentityIdentityx*
T0*'
_output_shapes
:€€€€€€€€€[

Identity_1IdentityIdentity:output:0*
T0*'
_output_shapes
:€€€€€€€€€§
	IdentityN	IdentityNIdentity:output:0x*
T
2*-
_gradient_op_typeCustomGradient-5671370*:
_output_shapes(
&:€€€€€€€€€:€€€€€€€€€\

Identity_2IdentityIdentityN:output:0*
T0*'
_output_shapes
:€€€€€€€€€"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*&
_input_shapes
:€€€€€€€€€:J F
'
_output_shapes
:€€€€€€€€€

_user_specified_namex
Ы

х
D__inference_dense_7_layer_call_and_return_conditional_losses_5670430

inputs0
matmul_readvariableop_resource:22-
biasadd_readvariableop_resource:2
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:22*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
«
Ц
)__inference_dense_6_layer_call_fn_5671418

inputs
unknown:22
	unknown_0:2
identityИҐStatefulPartitionedCallё
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_6_layer_call_and_return_conditional_losses_5670413o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
є
P
$__inference__update_step_xla_5671237
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
Ы

х
D__inference_dense_6_layer_call_and_return_conditional_losses_5670413

inputs0
matmul_readvariableop_resource:22-
biasadd_readvariableop_resource:2
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:22*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
Ѓ
л
)__inference_model_1_layer_call_fn_5670530
input_1
input_2
unknown:	2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
	unknown_5:2
	unknown_6:2
	unknown_7:22
	unknown_8:2
	unknown_9:22

unknown_10:2

unknown_11:22

unknown_12:2

unknown_13:22

unknown_14:2

unknown_15:2

unknown_16:
identity

identity_1ИҐStatefulPartitionedCall‘
StatefulPartitionedCallStatefulPartitionedCallinput_1input_2unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:€€€€€€€€€:€€€€€€€€€*4
_read_only_resource_inputs
	
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_model_1_layer_call_and_return_conditional_losses_5670489o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€q

Identity_1Identity StatefulPartitionedCall:output:1^NoOp*
T0*'
_output_shapes
:€€€€€€€€€`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*]
_input_shapesL
J:€€€€€€€€€	:€€€€€€€€€: : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:€€€€€€€€€	
!
_user_specified_name	input_1:PL
'
_output_shapes
:€€€€€€€€€
!
_user_specified_name	input_2
«
Ц
)__inference_dense_7_layer_call_fn_5671438

inputs
unknown:22
	unknown_0:2
identityИҐStatefulPartitionedCallё
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_7_layer_call_and_return_conditional_losses_5670430o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
ѓ
Y
-__inference_concatenate_layer_call_fn_5671382
inputs_0
inputs_1
identity≈
PartitionedCallPartitionedCallinputs_0inputs_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2 *0J 8В *Q
fLRJ
H__inference_concatenate_layer_call_and_return_conditional_losses_5670383`
IdentityIdentityPartitionedCall:output:0*
T0*'
_output_shapes
:€€€€€€€€€"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*9
_input_shapes(
&:€€€€€€€€€:€€€€€€€€€:Q M
'
_output_shapes
:€€€€€€€€€
"
_user_specified_name
inputs_0:QM
'
_output_shapes
:€€€€€€€€€
"
_user_specified_name
inputs_1
≠
L
$__inference__update_step_xla_5671262
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
Ѓ
л
)__inference_model_1_layer_call_fn_5670818
input_1
input_2
unknown:	2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
	unknown_5:2
	unknown_6:2
	unknown_7:22
	unknown_8:2
	unknown_9:22

unknown_10:2

unknown_11:22

unknown_12:2

unknown_13:22

unknown_14:2

unknown_15:2

unknown_16:
identity

identity_1ИҐStatefulPartitionedCall‘
StatefulPartitionedCallStatefulPartitionedCallinput_1input_2unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:€€€€€€€€€:€€€€€€€€€*4
_read_only_resource_inputs
	
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_model_1_layer_call_and_return_conditional_losses_5670733o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€q

Identity_1Identity StatefulPartitionedCall:output:1^NoOp*
T0*'
_output_shapes
:€€€€€€€€€`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*]
_input_shapesL
J:€€€€€€€€€	:€€€€€€€€€: : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:€€€€€€€€€	
!
_user_specified_name	input_1:PL
'
_output_shapes
:€€€€€€€€€
!
_user_specified_name	input_2
є
P
$__inference__update_step_xla_5671267
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
≠
L
$__inference__update_step_xla_5671252
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
Э
E
.__inference_grad_reverse_layer_call_fn_5671367
x
identityі
PartitionedCallPartitionedCallx*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2 *0J 8В *R
fMRK
I__inference_grad_reverse_layer_call_and_return_conditional_losses_5670374`
IdentityIdentityPartitionedCall:output:0*
T0*'
_output_shapes
:€€€€€€€€€"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*&
_input_shapes
:€€€€€€€€€:J F
'
_output_shapes
:€€€€€€€€€

_user_specified_namex
§R
•
D__inference_model_1_layer_call_and_return_conditional_losses_5671212
inputs_0
inputs_18
&dense_3_matmul_readvariableop_resource:	25
'dense_3_biasadd_readvariableop_resource:28
&dense_4_matmul_readvariableop_resource:225
'dense_4_biasadd_readvariableop_resource:24
"clf_matmul_readvariableop_resource:21
#clf_biasadd_readvariableop_resource:8
&dense_5_matmul_readvariableop_resource:25
'dense_5_biasadd_readvariableop_resource:28
&dense_6_matmul_readvariableop_resource:225
'dense_6_biasadd_readvariableop_resource:28
&dense_7_matmul_readvariableop_resource:225
'dense_7_biasadd_readvariableop_resource:28
&dense_8_matmul_readvariableop_resource:225
'dense_8_biasadd_readvariableop_resource:28
&dense_9_matmul_readvariableop_resource:225
'dense_9_biasadd_readvariableop_resource:24
"adv_matmul_readvariableop_resource:21
#adv_biasadd_readvariableop_resource:
identity

identity_1ИҐAdv/BiasAdd/ReadVariableOpҐAdv/MatMul/ReadVariableOpҐClf/BiasAdd/ReadVariableOpҐClf/MatMul/ReadVariableOpҐdense_3/BiasAdd/ReadVariableOpҐdense_3/MatMul/ReadVariableOpҐdense_4/BiasAdd/ReadVariableOpҐdense_4/MatMul/ReadVariableOpҐdense_5/BiasAdd/ReadVariableOpҐdense_5/MatMul/ReadVariableOpҐdense_6/BiasAdd/ReadVariableOpҐdense_6/MatMul/ReadVariableOpҐdense_7/BiasAdd/ReadVariableOpҐdense_7/MatMul/ReadVariableOpҐdense_8/BiasAdd/ReadVariableOpҐdense_8/MatMul/ReadVariableOpҐdense_9/BiasAdd/ReadVariableOpҐdense_9/MatMul/ReadVariableOpД
dense_3/MatMul/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource*
_output_shapes

:	2*
dtype0{
dense_3/MatMulMatMulinputs_0%dense_3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2В
dense_3/BiasAdd/ReadVariableOpReadVariableOp'dense_3_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0О
dense_3/BiasAddBiasAdddense_3/MatMul:product:0&dense_3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2`
dense_3/ReluReludense_3/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2Д
dense_4/MatMul/ReadVariableOpReadVariableOp&dense_4_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Н
dense_4/MatMulMatMuldense_3/Relu:activations:0%dense_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2В
dense_4/BiasAdd/ReadVariableOpReadVariableOp'dense_4_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0О
dense_4/BiasAddBiasAdddense_4/MatMul:product:0&dense_4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2`
dense_4/ReluReludense_4/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2|
Clf/MatMul/ReadVariableOpReadVariableOp"clf_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0Е

Clf/MatMulMatMuldense_4/Relu:activations:0!Clf/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€z
Clf/BiasAdd/ReadVariableOpReadVariableOp#clf_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0В
Clf/BiasAddBiasAddClf/MatMul:product:0"Clf/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€^
Clf/SigmoidSigmoidClf/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€d
grad_reverse/IdentityIdentityClf/Sigmoid:y:0*
T0*'
_output_shapes
:€€€€€€€€€u
grad_reverse/Identity_1Identitygrad_reverse/Identity:output:0*
T0*'
_output_shapes
:€€€€€€€€€ћ
grad_reverse/IdentityN	IdentityNgrad_reverse/Identity:output:0Clf/Sigmoid:y:0*
T
2*-
_gradient_op_typeCustomGradient-5671161*:
_output_shapes(
&:€€€€€€€€€:€€€€€€€€€Y
concatenate/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :¶
concatenate/concatConcatV2grad_reverse/IdentityN:output:0inputs_1 concatenate/concat/axis:output:0*
N*
T0*'
_output_shapes
:€€€€€€€€€Д
dense_5/MatMul/ReadVariableOpReadVariableOp&dense_5_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0О
dense_5/MatMulMatMulconcatenate/concat:output:0%dense_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2В
dense_5/BiasAdd/ReadVariableOpReadVariableOp'dense_5_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0О
dense_5/BiasAddBiasAdddense_5/MatMul:product:0&dense_5/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2`
dense_5/ReluReludense_5/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2Д
dense_6/MatMul/ReadVariableOpReadVariableOp&dense_6_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Н
dense_6/MatMulMatMuldense_5/Relu:activations:0%dense_6/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2В
dense_6/BiasAdd/ReadVariableOpReadVariableOp'dense_6_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0О
dense_6/BiasAddBiasAdddense_6/MatMul:product:0&dense_6/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2`
dense_6/ReluReludense_6/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2Д
dense_7/MatMul/ReadVariableOpReadVariableOp&dense_7_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Н
dense_7/MatMulMatMuldense_6/Relu:activations:0%dense_7/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2В
dense_7/BiasAdd/ReadVariableOpReadVariableOp'dense_7_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0О
dense_7/BiasAddBiasAdddense_7/MatMul:product:0&dense_7/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2`
dense_7/ReluReludense_7/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2Д
dense_8/MatMul/ReadVariableOpReadVariableOp&dense_8_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Н
dense_8/MatMulMatMuldense_7/Relu:activations:0%dense_8/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2В
dense_8/BiasAdd/ReadVariableOpReadVariableOp'dense_8_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0О
dense_8/BiasAddBiasAdddense_8/MatMul:product:0&dense_8/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2`
dense_8/ReluReludense_8/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2Д
dense_9/MatMul/ReadVariableOpReadVariableOp&dense_9_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Н
dense_9/MatMulMatMuldense_8/Relu:activations:0%dense_9/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2В
dense_9/BiasAdd/ReadVariableOpReadVariableOp'dense_9_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0О
dense_9/BiasAddBiasAdddense_9/MatMul:product:0&dense_9/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2`
dense_9/ReluReludense_9/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2|
Adv/MatMul/ReadVariableOpReadVariableOp"adv_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0Е

Adv/MatMulMatMuldense_9/Relu:activations:0!Adv/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€z
Adv/BiasAdd/ReadVariableOpReadVariableOp#adv_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0В
Adv/BiasAddBiasAddAdv/MatMul:product:0"Adv/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€^
Adv/SigmoidSigmoidAdv/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€^
IdentityIdentityClf/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€`

Identity_1IdentityAdv/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€€
NoOpNoOp^Adv/BiasAdd/ReadVariableOp^Adv/MatMul/ReadVariableOp^Clf/BiasAdd/ReadVariableOp^Clf/MatMul/ReadVariableOp^dense_3/BiasAdd/ReadVariableOp^dense_3/MatMul/ReadVariableOp^dense_4/BiasAdd/ReadVariableOp^dense_4/MatMul/ReadVariableOp^dense_5/BiasAdd/ReadVariableOp^dense_5/MatMul/ReadVariableOp^dense_6/BiasAdd/ReadVariableOp^dense_6/MatMul/ReadVariableOp^dense_7/BiasAdd/ReadVariableOp^dense_7/MatMul/ReadVariableOp^dense_8/BiasAdd/ReadVariableOp^dense_8/MatMul/ReadVariableOp^dense_9/BiasAdd/ReadVariableOp^dense_9/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*]
_input_shapesL
J:€€€€€€€€€	:€€€€€€€€€: : : : : : : : : : : : : : : : : : 28
Adv/BiasAdd/ReadVariableOpAdv/BiasAdd/ReadVariableOp26
Adv/MatMul/ReadVariableOpAdv/MatMul/ReadVariableOp28
Clf/BiasAdd/ReadVariableOpClf/BiasAdd/ReadVariableOp26
Clf/MatMul/ReadVariableOpClf/MatMul/ReadVariableOp2@
dense_3/BiasAdd/ReadVariableOpdense_3/BiasAdd/ReadVariableOp2>
dense_3/MatMul/ReadVariableOpdense_3/MatMul/ReadVariableOp2@
dense_4/BiasAdd/ReadVariableOpdense_4/BiasAdd/ReadVariableOp2>
dense_4/MatMul/ReadVariableOpdense_4/MatMul/ReadVariableOp2@
dense_5/BiasAdd/ReadVariableOpdense_5/BiasAdd/ReadVariableOp2>
dense_5/MatMul/ReadVariableOpdense_5/MatMul/ReadVariableOp2@
dense_6/BiasAdd/ReadVariableOpdense_6/BiasAdd/ReadVariableOp2>
dense_6/MatMul/ReadVariableOpdense_6/MatMul/ReadVariableOp2@
dense_7/BiasAdd/ReadVariableOpdense_7/BiasAdd/ReadVariableOp2>
dense_7/MatMul/ReadVariableOpdense_7/MatMul/ReadVariableOp2@
dense_8/BiasAdd/ReadVariableOpdense_8/BiasAdd/ReadVariableOp2>
dense_8/MatMul/ReadVariableOpdense_8/MatMul/ReadVariableOp2@
dense_9/BiasAdd/ReadVariableOpdense_9/BiasAdd/ReadVariableOp2>
dense_9/MatMul/ReadVariableOpdense_9/MatMul/ReadVariableOp:Q M
'
_output_shapes
:€€€€€€€€€	
"
_user_specified_name
inputs_0:QM
'
_output_shapes
:€€€€€€€€€
"
_user_specified_name
inputs_1
Ы

х
D__inference_dense_7_layer_call_and_return_conditional_losses_5671449

inputs0
matmul_readvariableop_resource:22-
biasadd_readvariableop_resource:2
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:22*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
є
P
$__inference__update_step_xla_5671247
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
Ы

х
D__inference_dense_4_layer_call_and_return_conditional_losses_5671342

inputs0
matmul_readvariableop_resource:22-
biasadd_readvariableop_resource:2
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:22*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
«
Ц
)__inference_dense_9_layer_call_fn_5671478

inputs
unknown:22
	unknown_0:2
identityИҐStatefulPartitionedCallё
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_9_layer_call_and_return_conditional_losses_5670464o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
«
Ц
)__inference_dense_4_layer_call_fn_5671331

inputs
unknown:22
	unknown_0:2
identityИҐStatefulPartitionedCallё
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_4_layer_call_and_return_conditional_losses_5670342o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
Ц

с
@__inference_Adv_layer_call_and_return_conditional_losses_5671509

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:2*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€V
SigmoidSigmoidBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€Z
IdentityIdentitySigmoid:y:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
Ы

х
D__inference_dense_3_layer_call_and_return_conditional_losses_5671322

inputs0
matmul_readvariableop_resource:	2-
biasadd_readvariableop_resource:2
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:	2*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€	: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€	
 
_user_specified_nameinputs
ђ
\
$__inference_internal_grad_fn_5671656
result_grads_0
result_grads_1
identityL
NegNegresult_grads_0*
T0*'
_output_shapes
:€€€€€€€€€J
mul/yConst*
_output_shapes
: *
dtype0*
valueB
 *   AU
mulMulNeg:y:0mul/y:output:0*
T0*'
_output_shapes
:€€€€€€€€€O
IdentityIdentitymul:z:0*
T0*'
_output_shapes
:€€€€€€€€€"
identityIdentity:output:0*9
_input_shapes(
&:€€€€€€€€€:€€€€€€€€€:W S
'
_output_shapes
:€€€€€€€€€
(
_user_specified_nameresult_grads_0:WS
'
_output_shapes
:€€€€€€€€€
(
_user_specified_nameresult_grads_1
њ
Т
%__inference_Clf_layer_call_fn_5671351

inputs
unknown:2
	unknown_0:
identityИҐStatefulPartitionedCallЏ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *I
fDRB
@__inference_Clf_layer_call_and_return_conditional_losses_5670359o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
≠
L
$__inference__update_step_xla_5671272
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
И
з
%__inference_signature_wrapper_5670972
input_1
input_2
unknown:	2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
	unknown_5:2
	unknown_6:2
	unknown_7:22
	unknown_8:2
	unknown_9:22

unknown_10:2

unknown_11:22

unknown_12:2

unknown_13:22

unknown_14:2

unknown_15:2

unknown_16:
identity

identity_1ИҐStatefulPartitionedCall≤
StatefulPartitionedCallStatefulPartitionedCallinput_1input_2unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:€€€€€€€€€:€€€€€€€€€*4
_read_only_resource_inputs
	
*2
config_proto" 

CPU

GPU2 *0J 8В *+
f&R$
"__inference__wrapped_model_5670305o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€q

Identity_1Identity StatefulPartitionedCall:output:1^NoOp*
T0*'
_output_shapes
:€€€€€€€€€`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*]
_input_shapesL
J:€€€€€€€€€	:€€€€€€€€€: : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:€€€€€€€€€	
!
_user_specified_name	input_1:PL
'
_output_shapes
:€€€€€€€€€
!
_user_specified_name	input_2
є
P
$__inference__update_step_xla_5671287
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
ђ
\
$__inference_internal_grad_fn_5671683
result_grads_0
result_grads_1
identityL
NegNegresult_grads_0*
T0*'
_output_shapes
:€€€€€€€€€J
mul/yConst*
_output_shapes
: *
dtype0*
valueB
 *   AU
mulMulNeg:y:0mul/y:output:0*
T0*'
_output_shapes
:€€€€€€€€€O
IdentityIdentitymul:z:0*
T0*'
_output_shapes
:€€€€€€€€€"
identityIdentity:output:0*9
_input_shapes(
&:€€€€€€€€€:€€€€€€€€€:W S
'
_output_shapes
:€€€€€€€€€
(
_user_specified_nameresult_grads_0:WS
'
_output_shapes
:€€€€€€€€€
(
_user_specified_nameresult_grads_1
є
P
$__inference__update_step_xla_5671217
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
«
Ц
)__inference_dense_5_layer_call_fn_5671398

inputs
unknown:2
	unknown_0:2
identityИҐStatefulPartitionedCallё
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_5_layer_call_and_return_conditional_losses_5670396o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:€€€€€€€€€
 
_user_specified_nameinputs
Ы

х
D__inference_dense_9_layer_call_and_return_conditional_losses_5671489

inputs0
matmul_readvariableop_resource:22-
biasadd_readvariableop_resource:2
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:22*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
≠
L
$__inference__update_step_xla_5671222
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
Е6
н
D__inference_model_1_layer_call_and_return_conditional_losses_5670924
input_1
input_2!
dense_3_5670875:	2
dense_3_5670877:2!
dense_4_5670880:22
dense_4_5670882:2
clf_5670885:2
clf_5670887:!
dense_5_5670892:2
dense_5_5670894:2!
dense_6_5670897:22
dense_6_5670899:2!
dense_7_5670902:22
dense_7_5670904:2!
dense_8_5670907:22
dense_8_5670909:2!
dense_9_5670912:22
dense_9_5670914:2
adv_5670917:2
adv_5670919:
identity

identity_1ИҐAdv/StatefulPartitionedCallҐClf/StatefulPartitionedCallҐdense_3/StatefulPartitionedCallҐdense_4/StatefulPartitionedCallҐdense_5/StatefulPartitionedCallҐdense_6/StatefulPartitionedCallҐdense_7/StatefulPartitionedCallҐdense_8/StatefulPartitionedCallҐdense_9/StatefulPartitionedCallх
dense_3/StatefulPartitionedCallStatefulPartitionedCallinput_1dense_3_5670875dense_3_5670877*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_3_layer_call_and_return_conditional_losses_5670325Ц
dense_4/StatefulPartitionedCallStatefulPartitionedCall(dense_3/StatefulPartitionedCall:output:0dense_4_5670880dense_4_5670882*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_4_layer_call_and_return_conditional_losses_5670342Ж
Clf/StatefulPartitionedCallStatefulPartitionedCall(dense_4/StatefulPartitionedCall:output:0clf_5670885clf_5670887*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *I
fDRB
@__inference_Clf_layer_call_and_return_conditional_losses_5670359д
grad_reverse/PartitionedCallPartitionedCall$Clf/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2 *0J 8В *R
fMRK
I__inference_grad_reverse_layer_call_and_return_conditional_losses_5670374н
concatenate/PartitionedCallPartitionedCall%grad_reverse/PartitionedCall:output:0input_2*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2 *0J 8В *Q
fLRJ
H__inference_concatenate_layer_call_and_return_conditional_losses_5670383Т
dense_5/StatefulPartitionedCallStatefulPartitionedCall$concatenate/PartitionedCall:output:0dense_5_5670892dense_5_5670894*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_5_layer_call_and_return_conditional_losses_5670396Ц
dense_6/StatefulPartitionedCallStatefulPartitionedCall(dense_5/StatefulPartitionedCall:output:0dense_6_5670897dense_6_5670899*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_6_layer_call_and_return_conditional_losses_5670413Ц
dense_7/StatefulPartitionedCallStatefulPartitionedCall(dense_6/StatefulPartitionedCall:output:0dense_7_5670902dense_7_5670904*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_7_layer_call_and_return_conditional_losses_5670430Ц
dense_8/StatefulPartitionedCallStatefulPartitionedCall(dense_7/StatefulPartitionedCall:output:0dense_8_5670907dense_8_5670909*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_8_layer_call_and_return_conditional_losses_5670447Ц
dense_9/StatefulPartitionedCallStatefulPartitionedCall(dense_8/StatefulPartitionedCall:output:0dense_9_5670912dense_9_5670914*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_9_layer_call_and_return_conditional_losses_5670464Ж
Adv/StatefulPartitionedCallStatefulPartitionedCall(dense_9/StatefulPartitionedCall:output:0adv_5670917adv_5670919*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *I
fDRB
@__inference_Adv_layer_call_and_return_conditional_losses_5670481s
IdentityIdentity$Clf/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€u

Identity_1Identity$Adv/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€р
NoOpNoOp^Adv/StatefulPartitionedCall^Clf/StatefulPartitionedCall ^dense_3/StatefulPartitionedCall ^dense_4/StatefulPartitionedCall ^dense_5/StatefulPartitionedCall ^dense_6/StatefulPartitionedCall ^dense_7/StatefulPartitionedCall ^dense_8/StatefulPartitionedCall ^dense_9/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*]
_input_shapesL
J:€€€€€€€€€	:€€€€€€€€€: : : : : : : : : : : : : : : : : : 2:
Adv/StatefulPartitionedCallAdv/StatefulPartitionedCall2:
Clf/StatefulPartitionedCallClf/StatefulPartitionedCall2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2B
dense_4/StatefulPartitionedCalldense_4/StatefulPartitionedCall2B
dense_5/StatefulPartitionedCalldense_5/StatefulPartitionedCall2B
dense_6/StatefulPartitionedCalldense_6/StatefulPartitionedCall2B
dense_7/StatefulPartitionedCalldense_7/StatefulPartitionedCall2B
dense_8/StatefulPartitionedCalldense_8/StatefulPartitionedCall2B
dense_9/StatefulPartitionedCalldense_9/StatefulPartitionedCall:P L
'
_output_shapes
:€€€€€€€€€	
!
_user_specified_name	input_1:PL
'
_output_shapes
:€€€€€€€€€
!
_user_specified_name	input_2
Г6
н
D__inference_model_1_layer_call_and_return_conditional_losses_5670489

inputs
inputs_1!
dense_3_5670326:	2
dense_3_5670328:2!
dense_4_5670343:22
dense_4_5670345:2
clf_5670360:2
clf_5670362:!
dense_5_5670397:2
dense_5_5670399:2!
dense_6_5670414:22
dense_6_5670416:2!
dense_7_5670431:22
dense_7_5670433:2!
dense_8_5670448:22
dense_8_5670450:2!
dense_9_5670465:22
dense_9_5670467:2
adv_5670482:2
adv_5670484:
identity

identity_1ИҐAdv/StatefulPartitionedCallҐClf/StatefulPartitionedCallҐdense_3/StatefulPartitionedCallҐdense_4/StatefulPartitionedCallҐdense_5/StatefulPartitionedCallҐdense_6/StatefulPartitionedCallҐdense_7/StatefulPartitionedCallҐdense_8/StatefulPartitionedCallҐdense_9/StatefulPartitionedCallф
dense_3/StatefulPartitionedCallStatefulPartitionedCallinputsdense_3_5670326dense_3_5670328*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_3_layer_call_and_return_conditional_losses_5670325Ц
dense_4/StatefulPartitionedCallStatefulPartitionedCall(dense_3/StatefulPartitionedCall:output:0dense_4_5670343dense_4_5670345*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_4_layer_call_and_return_conditional_losses_5670342Ж
Clf/StatefulPartitionedCallStatefulPartitionedCall(dense_4/StatefulPartitionedCall:output:0clf_5670360clf_5670362*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *I
fDRB
@__inference_Clf_layer_call_and_return_conditional_losses_5670359д
grad_reverse/PartitionedCallPartitionedCall$Clf/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2 *0J 8В *R
fMRK
I__inference_grad_reverse_layer_call_and_return_conditional_losses_5670374о
concatenate/PartitionedCallPartitionedCall%grad_reverse/PartitionedCall:output:0inputs_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2 *0J 8В *Q
fLRJ
H__inference_concatenate_layer_call_and_return_conditional_losses_5670383Т
dense_5/StatefulPartitionedCallStatefulPartitionedCall$concatenate/PartitionedCall:output:0dense_5_5670397dense_5_5670399*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_5_layer_call_and_return_conditional_losses_5670396Ц
dense_6/StatefulPartitionedCallStatefulPartitionedCall(dense_5/StatefulPartitionedCall:output:0dense_6_5670414dense_6_5670416*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_6_layer_call_and_return_conditional_losses_5670413Ц
dense_7/StatefulPartitionedCallStatefulPartitionedCall(dense_6/StatefulPartitionedCall:output:0dense_7_5670431dense_7_5670433*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_7_layer_call_and_return_conditional_losses_5670430Ц
dense_8/StatefulPartitionedCallStatefulPartitionedCall(dense_7/StatefulPartitionedCall:output:0dense_8_5670448dense_8_5670450*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_8_layer_call_and_return_conditional_losses_5670447Ц
dense_9/StatefulPartitionedCallStatefulPartitionedCall(dense_8/StatefulPartitionedCall:output:0dense_9_5670465dense_9_5670467*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_9_layer_call_and_return_conditional_losses_5670464Ж
Adv/StatefulPartitionedCallStatefulPartitionedCall(dense_9/StatefulPartitionedCall:output:0adv_5670482adv_5670484*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *I
fDRB
@__inference_Adv_layer_call_and_return_conditional_losses_5670481s
IdentityIdentity$Clf/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€u

Identity_1Identity$Adv/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€р
NoOpNoOp^Adv/StatefulPartitionedCall^Clf/StatefulPartitionedCall ^dense_3/StatefulPartitionedCall ^dense_4/StatefulPartitionedCall ^dense_5/StatefulPartitionedCall ^dense_6/StatefulPartitionedCall ^dense_7/StatefulPartitionedCall ^dense_8/StatefulPartitionedCall ^dense_9/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*]
_input_shapesL
J:€€€€€€€€€	:€€€€€€€€€: : : : : : : : : : : : : : : : : : 2:
Adv/StatefulPartitionedCallAdv/StatefulPartitionedCall2:
Clf/StatefulPartitionedCallClf/StatefulPartitionedCall2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2B
dense_4/StatefulPartitionedCalldense_4/StatefulPartitionedCall2B
dense_5/StatefulPartitionedCalldense_5/StatefulPartitionedCall2B
dense_6/StatefulPartitionedCalldense_6/StatefulPartitionedCall2B
dense_7/StatefulPartitionedCalldense_7/StatefulPartitionedCall2B
dense_8/StatefulPartitionedCalldense_8/StatefulPartitionedCall2B
dense_9/StatefulPartitionedCalldense_9/StatefulPartitionedCall:O K
'
_output_shapes
:€€€€€€€€€	
 
_user_specified_nameinputs:OK
'
_output_shapes
:€€€€€€€€€
 
_user_specified_nameinputs
§R
•
D__inference_model_1_layer_call_and_return_conditional_losses_5671136
inputs_0
inputs_18
&dense_3_matmul_readvariableop_resource:	25
'dense_3_biasadd_readvariableop_resource:28
&dense_4_matmul_readvariableop_resource:225
'dense_4_biasadd_readvariableop_resource:24
"clf_matmul_readvariableop_resource:21
#clf_biasadd_readvariableop_resource:8
&dense_5_matmul_readvariableop_resource:25
'dense_5_biasadd_readvariableop_resource:28
&dense_6_matmul_readvariableop_resource:225
'dense_6_biasadd_readvariableop_resource:28
&dense_7_matmul_readvariableop_resource:225
'dense_7_biasadd_readvariableop_resource:28
&dense_8_matmul_readvariableop_resource:225
'dense_8_biasadd_readvariableop_resource:28
&dense_9_matmul_readvariableop_resource:225
'dense_9_biasadd_readvariableop_resource:24
"adv_matmul_readvariableop_resource:21
#adv_biasadd_readvariableop_resource:
identity

identity_1ИҐAdv/BiasAdd/ReadVariableOpҐAdv/MatMul/ReadVariableOpҐClf/BiasAdd/ReadVariableOpҐClf/MatMul/ReadVariableOpҐdense_3/BiasAdd/ReadVariableOpҐdense_3/MatMul/ReadVariableOpҐdense_4/BiasAdd/ReadVariableOpҐdense_4/MatMul/ReadVariableOpҐdense_5/BiasAdd/ReadVariableOpҐdense_5/MatMul/ReadVariableOpҐdense_6/BiasAdd/ReadVariableOpҐdense_6/MatMul/ReadVariableOpҐdense_7/BiasAdd/ReadVariableOpҐdense_7/MatMul/ReadVariableOpҐdense_8/BiasAdd/ReadVariableOpҐdense_8/MatMul/ReadVariableOpҐdense_9/BiasAdd/ReadVariableOpҐdense_9/MatMul/ReadVariableOpД
dense_3/MatMul/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource*
_output_shapes

:	2*
dtype0{
dense_3/MatMulMatMulinputs_0%dense_3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2В
dense_3/BiasAdd/ReadVariableOpReadVariableOp'dense_3_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0О
dense_3/BiasAddBiasAdddense_3/MatMul:product:0&dense_3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2`
dense_3/ReluReludense_3/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2Д
dense_4/MatMul/ReadVariableOpReadVariableOp&dense_4_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Н
dense_4/MatMulMatMuldense_3/Relu:activations:0%dense_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2В
dense_4/BiasAdd/ReadVariableOpReadVariableOp'dense_4_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0О
dense_4/BiasAddBiasAdddense_4/MatMul:product:0&dense_4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2`
dense_4/ReluReludense_4/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2|
Clf/MatMul/ReadVariableOpReadVariableOp"clf_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0Е

Clf/MatMulMatMuldense_4/Relu:activations:0!Clf/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€z
Clf/BiasAdd/ReadVariableOpReadVariableOp#clf_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0В
Clf/BiasAddBiasAddClf/MatMul:product:0"Clf/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€^
Clf/SigmoidSigmoidClf/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€d
grad_reverse/IdentityIdentityClf/Sigmoid:y:0*
T0*'
_output_shapes
:€€€€€€€€€u
grad_reverse/Identity_1Identitygrad_reverse/Identity:output:0*
T0*'
_output_shapes
:€€€€€€€€€ћ
grad_reverse/IdentityN	IdentityNgrad_reverse/Identity:output:0Clf/Sigmoid:y:0*
T
2*-
_gradient_op_typeCustomGradient-5671085*:
_output_shapes(
&:€€€€€€€€€:€€€€€€€€€Y
concatenate/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :¶
concatenate/concatConcatV2grad_reverse/IdentityN:output:0inputs_1 concatenate/concat/axis:output:0*
N*
T0*'
_output_shapes
:€€€€€€€€€Д
dense_5/MatMul/ReadVariableOpReadVariableOp&dense_5_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0О
dense_5/MatMulMatMulconcatenate/concat:output:0%dense_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2В
dense_5/BiasAdd/ReadVariableOpReadVariableOp'dense_5_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0О
dense_5/BiasAddBiasAdddense_5/MatMul:product:0&dense_5/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2`
dense_5/ReluReludense_5/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2Д
dense_6/MatMul/ReadVariableOpReadVariableOp&dense_6_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Н
dense_6/MatMulMatMuldense_5/Relu:activations:0%dense_6/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2В
dense_6/BiasAdd/ReadVariableOpReadVariableOp'dense_6_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0О
dense_6/BiasAddBiasAdddense_6/MatMul:product:0&dense_6/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2`
dense_6/ReluReludense_6/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2Д
dense_7/MatMul/ReadVariableOpReadVariableOp&dense_7_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Н
dense_7/MatMulMatMuldense_6/Relu:activations:0%dense_7/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2В
dense_7/BiasAdd/ReadVariableOpReadVariableOp'dense_7_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0О
dense_7/BiasAddBiasAdddense_7/MatMul:product:0&dense_7/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2`
dense_7/ReluReludense_7/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2Д
dense_8/MatMul/ReadVariableOpReadVariableOp&dense_8_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Н
dense_8/MatMulMatMuldense_7/Relu:activations:0%dense_8/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2В
dense_8/BiasAdd/ReadVariableOpReadVariableOp'dense_8_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0О
dense_8/BiasAddBiasAdddense_8/MatMul:product:0&dense_8/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2`
dense_8/ReluReludense_8/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2Д
dense_9/MatMul/ReadVariableOpReadVariableOp&dense_9_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0Н
dense_9/MatMulMatMuldense_8/Relu:activations:0%dense_9/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2В
dense_9/BiasAdd/ReadVariableOpReadVariableOp'dense_9_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0О
dense_9/BiasAddBiasAdddense_9/MatMul:product:0&dense_9/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2`
dense_9/ReluReludense_9/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2|
Adv/MatMul/ReadVariableOpReadVariableOp"adv_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0Е

Adv/MatMulMatMuldense_9/Relu:activations:0!Adv/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€z
Adv/BiasAdd/ReadVariableOpReadVariableOp#adv_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0В
Adv/BiasAddBiasAddAdv/MatMul:product:0"Adv/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€^
Adv/SigmoidSigmoidAdv/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€^
IdentityIdentityClf/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€`

Identity_1IdentityAdv/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€€
NoOpNoOp^Adv/BiasAdd/ReadVariableOp^Adv/MatMul/ReadVariableOp^Clf/BiasAdd/ReadVariableOp^Clf/MatMul/ReadVariableOp^dense_3/BiasAdd/ReadVariableOp^dense_3/MatMul/ReadVariableOp^dense_4/BiasAdd/ReadVariableOp^dense_4/MatMul/ReadVariableOp^dense_5/BiasAdd/ReadVariableOp^dense_5/MatMul/ReadVariableOp^dense_6/BiasAdd/ReadVariableOp^dense_6/MatMul/ReadVariableOp^dense_7/BiasAdd/ReadVariableOp^dense_7/MatMul/ReadVariableOp^dense_8/BiasAdd/ReadVariableOp^dense_8/MatMul/ReadVariableOp^dense_9/BiasAdd/ReadVariableOp^dense_9/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*]
_input_shapesL
J:€€€€€€€€€	:€€€€€€€€€: : : : : : : : : : : : : : : : : : 28
Adv/BiasAdd/ReadVariableOpAdv/BiasAdd/ReadVariableOp26
Adv/MatMul/ReadVariableOpAdv/MatMul/ReadVariableOp28
Clf/BiasAdd/ReadVariableOpClf/BiasAdd/ReadVariableOp26
Clf/MatMul/ReadVariableOpClf/MatMul/ReadVariableOp2@
dense_3/BiasAdd/ReadVariableOpdense_3/BiasAdd/ReadVariableOp2>
dense_3/MatMul/ReadVariableOpdense_3/MatMul/ReadVariableOp2@
dense_4/BiasAdd/ReadVariableOpdense_4/BiasAdd/ReadVariableOp2>
dense_4/MatMul/ReadVariableOpdense_4/MatMul/ReadVariableOp2@
dense_5/BiasAdd/ReadVariableOpdense_5/BiasAdd/ReadVariableOp2>
dense_5/MatMul/ReadVariableOpdense_5/MatMul/ReadVariableOp2@
dense_6/BiasAdd/ReadVariableOpdense_6/BiasAdd/ReadVariableOp2>
dense_6/MatMul/ReadVariableOpdense_6/MatMul/ReadVariableOp2@
dense_7/BiasAdd/ReadVariableOpdense_7/BiasAdd/ReadVariableOp2>
dense_7/MatMul/ReadVariableOpdense_7/MatMul/ReadVariableOp2@
dense_8/BiasAdd/ReadVariableOpdense_8/BiasAdd/ReadVariableOp2>
dense_8/MatMul/ReadVariableOpdense_8/MatMul/ReadVariableOp2@
dense_9/BiasAdd/ReadVariableOpdense_9/BiasAdd/ReadVariableOp2>
dense_9/MatMul/ReadVariableOpdense_9/MatMul/ReadVariableOp:Q M
'
_output_shapes
:€€€€€€€€€	
"
_user_specified_name
inputs_0:QM
'
_output_shapes
:€€€€€€€€€
"
_user_specified_name
inputs_1
≠
L
$__inference__update_step_xla_5671232
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
Г6
н
D__inference_model_1_layer_call_and_return_conditional_losses_5670733

inputs
inputs_1!
dense_3_5670684:	2
dense_3_5670686:2!
dense_4_5670689:22
dense_4_5670691:2
clf_5670694:2
clf_5670696:!
dense_5_5670701:2
dense_5_5670703:2!
dense_6_5670706:22
dense_6_5670708:2!
dense_7_5670711:22
dense_7_5670713:2!
dense_8_5670716:22
dense_8_5670718:2!
dense_9_5670721:22
dense_9_5670723:2
adv_5670726:2
adv_5670728:
identity

identity_1ИҐAdv/StatefulPartitionedCallҐClf/StatefulPartitionedCallҐdense_3/StatefulPartitionedCallҐdense_4/StatefulPartitionedCallҐdense_5/StatefulPartitionedCallҐdense_6/StatefulPartitionedCallҐdense_7/StatefulPartitionedCallҐdense_8/StatefulPartitionedCallҐdense_9/StatefulPartitionedCallф
dense_3/StatefulPartitionedCallStatefulPartitionedCallinputsdense_3_5670684dense_3_5670686*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_3_layer_call_and_return_conditional_losses_5670325Ц
dense_4/StatefulPartitionedCallStatefulPartitionedCall(dense_3/StatefulPartitionedCall:output:0dense_4_5670689dense_4_5670691*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_4_layer_call_and_return_conditional_losses_5670342Ж
Clf/StatefulPartitionedCallStatefulPartitionedCall(dense_4/StatefulPartitionedCall:output:0clf_5670694clf_5670696*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *I
fDRB
@__inference_Clf_layer_call_and_return_conditional_losses_5670359д
grad_reverse/PartitionedCallPartitionedCall$Clf/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2 *0J 8В *R
fMRK
I__inference_grad_reverse_layer_call_and_return_conditional_losses_5670374о
concatenate/PartitionedCallPartitionedCall%grad_reverse/PartitionedCall:output:0inputs_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2 *0J 8В *Q
fLRJ
H__inference_concatenate_layer_call_and_return_conditional_losses_5670383Т
dense_5/StatefulPartitionedCallStatefulPartitionedCall$concatenate/PartitionedCall:output:0dense_5_5670701dense_5_5670703*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_5_layer_call_and_return_conditional_losses_5670396Ц
dense_6/StatefulPartitionedCallStatefulPartitionedCall(dense_5/StatefulPartitionedCall:output:0dense_6_5670706dense_6_5670708*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_6_layer_call_and_return_conditional_losses_5670413Ц
dense_7/StatefulPartitionedCallStatefulPartitionedCall(dense_6/StatefulPartitionedCall:output:0dense_7_5670711dense_7_5670713*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_7_layer_call_and_return_conditional_losses_5670430Ц
dense_8/StatefulPartitionedCallStatefulPartitionedCall(dense_7/StatefulPartitionedCall:output:0dense_8_5670716dense_8_5670718*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_8_layer_call_and_return_conditional_losses_5670447Ц
dense_9/StatefulPartitionedCallStatefulPartitionedCall(dense_8/StatefulPartitionedCall:output:0dense_9_5670721dense_9_5670723*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_9_layer_call_and_return_conditional_losses_5670464Ж
Adv/StatefulPartitionedCallStatefulPartitionedCall(dense_9/StatefulPartitionedCall:output:0adv_5670726adv_5670728*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *I
fDRB
@__inference_Adv_layer_call_and_return_conditional_losses_5670481s
IdentityIdentity$Clf/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€u

Identity_1Identity$Adv/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€р
NoOpNoOp^Adv/StatefulPartitionedCall^Clf/StatefulPartitionedCall ^dense_3/StatefulPartitionedCall ^dense_4/StatefulPartitionedCall ^dense_5/StatefulPartitionedCall ^dense_6/StatefulPartitionedCall ^dense_7/StatefulPartitionedCall ^dense_8/StatefulPartitionedCall ^dense_9/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*]
_input_shapesL
J:€€€€€€€€€	:€€€€€€€€€: : : : : : : : : : : : : : : : : : 2:
Adv/StatefulPartitionedCallAdv/StatefulPartitionedCall2:
Clf/StatefulPartitionedCallClf/StatefulPartitionedCall2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2B
dense_4/StatefulPartitionedCalldense_4/StatefulPartitionedCall2B
dense_5/StatefulPartitionedCalldense_5/StatefulPartitionedCall2B
dense_6/StatefulPartitionedCalldense_6/StatefulPartitionedCall2B
dense_7/StatefulPartitionedCalldense_7/StatefulPartitionedCall2B
dense_8/StatefulPartitionedCalldense_8/StatefulPartitionedCall2B
dense_9/StatefulPartitionedCalldense_9/StatefulPartitionedCall:O K
'
_output_shapes
:€€€€€€€€€	
 
_user_specified_nameinputs:OK
'
_output_shapes
:€€€€€€€€€
 
_user_specified_nameinputs
Ы

х
D__inference_dense_5_layer_call_and_return_conditional_losses_5671409

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:2
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:2*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€
 
_user_specified_nameinputs
Ы

х
D__inference_dense_3_layer_call_and_return_conditional_losses_5670325

inputs0
matmul_readvariableop_resource:	2-
biasadd_readvariableop_resource:2
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:	2*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€	: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€	
 
_user_specified_nameinputs
є
r
H__inference_concatenate_layer_call_and_return_conditional_losses_5670383

inputs
inputs_1
identityM
concat/axisConst*
_output_shapes
: *
dtype0*
value	B :u
concatConcatV2inputsinputs_1concat/axis:output:0*
N*
T0*'
_output_shapes
:€€€€€€€€€W
IdentityIdentityconcat:output:0*
T0*'
_output_shapes
:€€€€€€€€€"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*9
_input_shapes(
&:€€€€€€€€€:€€€€€€€€€:O K
'
_output_shapes
:€€€€€€€€€
 
_user_specified_nameinputs:OK
'
_output_shapes
:€€€€€€€€€
 
_user_specified_nameinputs
«
Ц
)__inference_dense_8_layer_call_fn_5671458

inputs
unknown:22
	unknown_0:2
identityИҐStatefulPartitionedCallё
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_8_layer_call_and_return_conditional_losses_5670447o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
і
н
)__inference_model_1_layer_call_fn_5671060
inputs_0
inputs_1
unknown:	2
	unknown_0:2
	unknown_1:22
	unknown_2:2
	unknown_3:2
	unknown_4:
	unknown_5:2
	unknown_6:2
	unknown_7:22
	unknown_8:2
	unknown_9:22

unknown_10:2

unknown_11:22

unknown_12:2

unknown_13:22

unknown_14:2

unknown_15:2

unknown_16:
identity

identity_1ИҐStatefulPartitionedCall÷
StatefulPartitionedCallStatefulPartitionedCallinputs_0inputs_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:€€€€€€€€€:€€€€€€€€€*4
_read_only_resource_inputs
	
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_model_1_layer_call_and_return_conditional_losses_5670733o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€q

Identity_1Identity StatefulPartitionedCall:output:1^NoOp*
T0*'
_output_shapes
:€€€€€€€€€`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*]
_input_shapesL
J:€€€€€€€€€	:€€€€€€€€€: : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:Q M
'
_output_shapes
:€€€€€€€€€	
"
_user_specified_name
inputs_0:QM
'
_output_shapes
:€€€€€€€€€
"
_user_specified_name
inputs_1
Т_
°
"__inference__wrapped_model_5670305
input_1
input_2@
.model_1_dense_3_matmul_readvariableop_resource:	2=
/model_1_dense_3_biasadd_readvariableop_resource:2@
.model_1_dense_4_matmul_readvariableop_resource:22=
/model_1_dense_4_biasadd_readvariableop_resource:2<
*model_1_clf_matmul_readvariableop_resource:29
+model_1_clf_biasadd_readvariableop_resource:@
.model_1_dense_5_matmul_readvariableop_resource:2=
/model_1_dense_5_biasadd_readvariableop_resource:2@
.model_1_dense_6_matmul_readvariableop_resource:22=
/model_1_dense_6_biasadd_readvariableop_resource:2@
.model_1_dense_7_matmul_readvariableop_resource:22=
/model_1_dense_7_biasadd_readvariableop_resource:2@
.model_1_dense_8_matmul_readvariableop_resource:22=
/model_1_dense_8_biasadd_readvariableop_resource:2@
.model_1_dense_9_matmul_readvariableop_resource:22=
/model_1_dense_9_biasadd_readvariableop_resource:2<
*model_1_adv_matmul_readvariableop_resource:29
+model_1_adv_biasadd_readvariableop_resource:
identity

identity_1ИҐ"model_1/Adv/BiasAdd/ReadVariableOpҐ!model_1/Adv/MatMul/ReadVariableOpҐ"model_1/Clf/BiasAdd/ReadVariableOpҐ!model_1/Clf/MatMul/ReadVariableOpҐ&model_1/dense_3/BiasAdd/ReadVariableOpҐ%model_1/dense_3/MatMul/ReadVariableOpҐ&model_1/dense_4/BiasAdd/ReadVariableOpҐ%model_1/dense_4/MatMul/ReadVariableOpҐ&model_1/dense_5/BiasAdd/ReadVariableOpҐ%model_1/dense_5/MatMul/ReadVariableOpҐ&model_1/dense_6/BiasAdd/ReadVariableOpҐ%model_1/dense_6/MatMul/ReadVariableOpҐ&model_1/dense_7/BiasAdd/ReadVariableOpҐ%model_1/dense_7/MatMul/ReadVariableOpҐ&model_1/dense_8/BiasAdd/ReadVariableOpҐ%model_1/dense_8/MatMul/ReadVariableOpҐ&model_1/dense_9/BiasAdd/ReadVariableOpҐ%model_1/dense_9/MatMul/ReadVariableOpФ
%model_1/dense_3/MatMul/ReadVariableOpReadVariableOp.model_1_dense_3_matmul_readvariableop_resource*
_output_shapes

:	2*
dtype0К
model_1/dense_3/MatMulMatMulinput_1-model_1/dense_3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2Т
&model_1/dense_3/BiasAdd/ReadVariableOpReadVariableOp/model_1_dense_3_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0¶
model_1/dense_3/BiasAddBiasAdd model_1/dense_3/MatMul:product:0.model_1/dense_3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2p
model_1/dense_3/ReluRelu model_1/dense_3/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2Ф
%model_1/dense_4/MatMul/ReadVariableOpReadVariableOp.model_1_dense_4_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0•
model_1/dense_4/MatMulMatMul"model_1/dense_3/Relu:activations:0-model_1/dense_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2Т
&model_1/dense_4/BiasAdd/ReadVariableOpReadVariableOp/model_1_dense_4_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0¶
model_1/dense_4/BiasAddBiasAdd model_1/dense_4/MatMul:product:0.model_1/dense_4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2p
model_1/dense_4/ReluRelu model_1/dense_4/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2М
!model_1/Clf/MatMul/ReadVariableOpReadVariableOp*model_1_clf_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0Э
model_1/Clf/MatMulMatMul"model_1/dense_4/Relu:activations:0)model_1/Clf/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€К
"model_1/Clf/BiasAdd/ReadVariableOpReadVariableOp+model_1_clf_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0Ъ
model_1/Clf/BiasAddBiasAddmodel_1/Clf/MatMul:product:0*model_1/Clf/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€n
model_1/Clf/SigmoidSigmoidmodel_1/Clf/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€t
model_1/grad_reverse/IdentityIdentitymodel_1/Clf/Sigmoid:y:0*
T0*'
_output_shapes
:€€€€€€€€€Е
model_1/grad_reverse/Identity_1Identity&model_1/grad_reverse/Identity:output:0*
T0*'
_output_shapes
:€€€€€€€€€д
model_1/grad_reverse/IdentityN	IdentityN&model_1/grad_reverse/Identity:output:0model_1/Clf/Sigmoid:y:0*
T
2*-
_gradient_op_typeCustomGradient-5670254*:
_output_shapes(
&:€€€€€€€€€:€€€€€€€€€a
model_1/concatenate/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :љ
model_1/concatenate/concatConcatV2'model_1/grad_reverse/IdentityN:output:0input_2(model_1/concatenate/concat/axis:output:0*
N*
T0*'
_output_shapes
:€€€€€€€€€Ф
%model_1/dense_5/MatMul/ReadVariableOpReadVariableOp.model_1_dense_5_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0¶
model_1/dense_5/MatMulMatMul#model_1/concatenate/concat:output:0-model_1/dense_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2Т
&model_1/dense_5/BiasAdd/ReadVariableOpReadVariableOp/model_1_dense_5_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0¶
model_1/dense_5/BiasAddBiasAdd model_1/dense_5/MatMul:product:0.model_1/dense_5/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2p
model_1/dense_5/ReluRelu model_1/dense_5/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2Ф
%model_1/dense_6/MatMul/ReadVariableOpReadVariableOp.model_1_dense_6_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0•
model_1/dense_6/MatMulMatMul"model_1/dense_5/Relu:activations:0-model_1/dense_6/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2Т
&model_1/dense_6/BiasAdd/ReadVariableOpReadVariableOp/model_1_dense_6_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0¶
model_1/dense_6/BiasAddBiasAdd model_1/dense_6/MatMul:product:0.model_1/dense_6/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2p
model_1/dense_6/ReluRelu model_1/dense_6/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2Ф
%model_1/dense_7/MatMul/ReadVariableOpReadVariableOp.model_1_dense_7_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0•
model_1/dense_7/MatMulMatMul"model_1/dense_6/Relu:activations:0-model_1/dense_7/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2Т
&model_1/dense_7/BiasAdd/ReadVariableOpReadVariableOp/model_1_dense_7_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0¶
model_1/dense_7/BiasAddBiasAdd model_1/dense_7/MatMul:product:0.model_1/dense_7/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2p
model_1/dense_7/ReluRelu model_1/dense_7/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2Ф
%model_1/dense_8/MatMul/ReadVariableOpReadVariableOp.model_1_dense_8_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0•
model_1/dense_8/MatMulMatMul"model_1/dense_7/Relu:activations:0-model_1/dense_8/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2Т
&model_1/dense_8/BiasAdd/ReadVariableOpReadVariableOp/model_1_dense_8_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0¶
model_1/dense_8/BiasAddBiasAdd model_1/dense_8/MatMul:product:0.model_1/dense_8/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2p
model_1/dense_8/ReluRelu model_1/dense_8/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2Ф
%model_1/dense_9/MatMul/ReadVariableOpReadVariableOp.model_1_dense_9_matmul_readvariableop_resource*
_output_shapes

:22*
dtype0•
model_1/dense_9/MatMulMatMul"model_1/dense_8/Relu:activations:0-model_1/dense_9/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2Т
&model_1/dense_9/BiasAdd/ReadVariableOpReadVariableOp/model_1_dense_9_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype0¶
model_1/dense_9/BiasAddBiasAdd model_1/dense_9/MatMul:product:0.model_1/dense_9/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2p
model_1/dense_9/ReluRelu model_1/dense_9/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2М
!model_1/Adv/MatMul/ReadVariableOpReadVariableOp*model_1_adv_matmul_readvariableop_resource*
_output_shapes

:2*
dtype0Э
model_1/Adv/MatMulMatMul"model_1/dense_9/Relu:activations:0)model_1/Adv/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€К
"model_1/Adv/BiasAdd/ReadVariableOpReadVariableOp+model_1_adv_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0Ъ
model_1/Adv/BiasAddBiasAddmodel_1/Adv/MatMul:product:0*model_1/Adv/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€n
model_1/Adv/SigmoidSigmoidmodel_1/Adv/BiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€f
IdentityIdentitymodel_1/Adv/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€h

Identity_1Identitymodel_1/Clf/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€П
NoOpNoOp#^model_1/Adv/BiasAdd/ReadVariableOp"^model_1/Adv/MatMul/ReadVariableOp#^model_1/Clf/BiasAdd/ReadVariableOp"^model_1/Clf/MatMul/ReadVariableOp'^model_1/dense_3/BiasAdd/ReadVariableOp&^model_1/dense_3/MatMul/ReadVariableOp'^model_1/dense_4/BiasAdd/ReadVariableOp&^model_1/dense_4/MatMul/ReadVariableOp'^model_1/dense_5/BiasAdd/ReadVariableOp&^model_1/dense_5/MatMul/ReadVariableOp'^model_1/dense_6/BiasAdd/ReadVariableOp&^model_1/dense_6/MatMul/ReadVariableOp'^model_1/dense_7/BiasAdd/ReadVariableOp&^model_1/dense_7/MatMul/ReadVariableOp'^model_1/dense_8/BiasAdd/ReadVariableOp&^model_1/dense_8/MatMul/ReadVariableOp'^model_1/dense_9/BiasAdd/ReadVariableOp&^model_1/dense_9/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*]
_input_shapesL
J:€€€€€€€€€	:€€€€€€€€€: : : : : : : : : : : : : : : : : : 2H
"model_1/Adv/BiasAdd/ReadVariableOp"model_1/Adv/BiasAdd/ReadVariableOp2F
!model_1/Adv/MatMul/ReadVariableOp!model_1/Adv/MatMul/ReadVariableOp2H
"model_1/Clf/BiasAdd/ReadVariableOp"model_1/Clf/BiasAdd/ReadVariableOp2F
!model_1/Clf/MatMul/ReadVariableOp!model_1/Clf/MatMul/ReadVariableOp2P
&model_1/dense_3/BiasAdd/ReadVariableOp&model_1/dense_3/BiasAdd/ReadVariableOp2N
%model_1/dense_3/MatMul/ReadVariableOp%model_1/dense_3/MatMul/ReadVariableOp2P
&model_1/dense_4/BiasAdd/ReadVariableOp&model_1/dense_4/BiasAdd/ReadVariableOp2N
%model_1/dense_4/MatMul/ReadVariableOp%model_1/dense_4/MatMul/ReadVariableOp2P
&model_1/dense_5/BiasAdd/ReadVariableOp&model_1/dense_5/BiasAdd/ReadVariableOp2N
%model_1/dense_5/MatMul/ReadVariableOp%model_1/dense_5/MatMul/ReadVariableOp2P
&model_1/dense_6/BiasAdd/ReadVariableOp&model_1/dense_6/BiasAdd/ReadVariableOp2N
%model_1/dense_6/MatMul/ReadVariableOp%model_1/dense_6/MatMul/ReadVariableOp2P
&model_1/dense_7/BiasAdd/ReadVariableOp&model_1/dense_7/BiasAdd/ReadVariableOp2N
%model_1/dense_7/MatMul/ReadVariableOp%model_1/dense_7/MatMul/ReadVariableOp2P
&model_1/dense_8/BiasAdd/ReadVariableOp&model_1/dense_8/BiasAdd/ReadVariableOp2N
%model_1/dense_8/MatMul/ReadVariableOp%model_1/dense_8/MatMul/ReadVariableOp2P
&model_1/dense_9/BiasAdd/ReadVariableOp&model_1/dense_9/BiasAdd/ReadVariableOp2N
%model_1/dense_9/MatMul/ReadVariableOp%model_1/dense_9/MatMul/ReadVariableOp:P L
'
_output_shapes
:€€€€€€€€€	
!
_user_specified_name	input_1:PL
'
_output_shapes
:€€€€€€€€€
!
_user_specified_name	input_2
≠
L
$__inference__update_step_xla_5671282
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
≠
L
$__inference__update_step_xla_5671302
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
ђ
\
$__inference_internal_grad_fn_5671674
result_grads_0
result_grads_1
identityL
NegNegresult_grads_0*
T0*'
_output_shapes
:€€€€€€€€€J
mul/yConst*
_output_shapes
: *
dtype0*
valueB
 *   AU
mulMulNeg:y:0mul/y:output:0*
T0*'
_output_shapes
:€€€€€€€€€O
IdentityIdentitymul:z:0*
T0*'
_output_shapes
:€€€€€€€€€"
identityIdentity:output:0*9
_input_shapes(
&:€€€€€€€€€:€€€€€€€€€:W S
'
_output_shapes
:€€€€€€€€€
(
_user_specified_nameresult_grads_0:WS
'
_output_shapes
:€€€€€€€€€
(
_user_specified_nameresult_grads_1
≠
L
$__inference__update_step_xla_5671242
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
є
P
$__inference__update_step_xla_5671297
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
Е6
н
D__inference_model_1_layer_call_and_return_conditional_losses_5670871
input_1
input_2!
dense_3_5670822:	2
dense_3_5670824:2!
dense_4_5670827:22
dense_4_5670829:2
clf_5670832:2
clf_5670834:!
dense_5_5670839:2
dense_5_5670841:2!
dense_6_5670844:22
dense_6_5670846:2!
dense_7_5670849:22
dense_7_5670851:2!
dense_8_5670854:22
dense_8_5670856:2!
dense_9_5670859:22
dense_9_5670861:2
adv_5670864:2
adv_5670866:
identity

identity_1ИҐAdv/StatefulPartitionedCallҐClf/StatefulPartitionedCallҐdense_3/StatefulPartitionedCallҐdense_4/StatefulPartitionedCallҐdense_5/StatefulPartitionedCallҐdense_6/StatefulPartitionedCallҐdense_7/StatefulPartitionedCallҐdense_8/StatefulPartitionedCallҐdense_9/StatefulPartitionedCallх
dense_3/StatefulPartitionedCallStatefulPartitionedCallinput_1dense_3_5670822dense_3_5670824*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_3_layer_call_and_return_conditional_losses_5670325Ц
dense_4/StatefulPartitionedCallStatefulPartitionedCall(dense_3/StatefulPartitionedCall:output:0dense_4_5670827dense_4_5670829*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_4_layer_call_and_return_conditional_losses_5670342Ж
Clf/StatefulPartitionedCallStatefulPartitionedCall(dense_4/StatefulPartitionedCall:output:0clf_5670832clf_5670834*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *I
fDRB
@__inference_Clf_layer_call_and_return_conditional_losses_5670359д
grad_reverse/PartitionedCallPartitionedCall$Clf/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2 *0J 8В *R
fMRK
I__inference_grad_reverse_layer_call_and_return_conditional_losses_5670374н
concatenate/PartitionedCallPartitionedCall%grad_reverse/PartitionedCall:output:0input_2*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2 *0J 8В *Q
fLRJ
H__inference_concatenate_layer_call_and_return_conditional_losses_5670383Т
dense_5/StatefulPartitionedCallStatefulPartitionedCall$concatenate/PartitionedCall:output:0dense_5_5670839dense_5_5670841*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_5_layer_call_and_return_conditional_losses_5670396Ц
dense_6/StatefulPartitionedCallStatefulPartitionedCall(dense_5/StatefulPartitionedCall:output:0dense_6_5670844dense_6_5670846*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_6_layer_call_and_return_conditional_losses_5670413Ц
dense_7/StatefulPartitionedCallStatefulPartitionedCall(dense_6/StatefulPartitionedCall:output:0dense_7_5670849dense_7_5670851*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_7_layer_call_and_return_conditional_losses_5670430Ц
dense_8/StatefulPartitionedCallStatefulPartitionedCall(dense_7/StatefulPartitionedCall:output:0dense_8_5670854dense_8_5670856*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_8_layer_call_and_return_conditional_losses_5670447Ц
dense_9/StatefulPartitionedCallStatefulPartitionedCall(dense_8/StatefulPartitionedCall:output:0dense_9_5670859dense_9_5670861*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€2*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *M
fHRF
D__inference_dense_9_layer_call_and_return_conditional_losses_5670464Ж
Adv/StatefulPartitionedCallStatefulPartitionedCall(dense_9/StatefulPartitionedCall:output:0adv_5670864adv_5670866*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:€€€€€€€€€*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2 *0J 8В *I
fDRB
@__inference_Adv_layer_call_and_return_conditional_losses_5670481s
IdentityIdentity$Clf/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€u

Identity_1Identity$Adv/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€р
NoOpNoOp^Adv/StatefulPartitionedCall^Clf/StatefulPartitionedCall ^dense_3/StatefulPartitionedCall ^dense_4/StatefulPartitionedCall ^dense_5/StatefulPartitionedCall ^dense_6/StatefulPartitionedCall ^dense_7/StatefulPartitionedCall ^dense_8/StatefulPartitionedCall ^dense_9/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*]
_input_shapesL
J:€€€€€€€€€	:€€€€€€€€€: : : : : : : : : : : : : : : : : : 2:
Adv/StatefulPartitionedCallAdv/StatefulPartitionedCall2:
Clf/StatefulPartitionedCallClf/StatefulPartitionedCall2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2B
dense_4/StatefulPartitionedCalldense_4/StatefulPartitionedCall2B
dense_5/StatefulPartitionedCalldense_5/StatefulPartitionedCall2B
dense_6/StatefulPartitionedCalldense_6/StatefulPartitionedCall2B
dense_7/StatefulPartitionedCalldense_7/StatefulPartitionedCall2B
dense_8/StatefulPartitionedCalldense_8/StatefulPartitionedCall2B
dense_9/StatefulPartitionedCalldense_9/StatefulPartitionedCall:P L
'
_output_shapes
:€€€€€€€€€	
!
_user_specified_name	input_1:PL
'
_output_shapes
:€€€€€€€€€
!
_user_specified_name	input_2
Ы

х
D__inference_dense_4_layer_call_and_return_conditional_losses_5670342

inputs0
matmul_readvariableop_resource:22-
biasadd_readvariableop_resource:2
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:22*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
Ц

с
@__inference_Adv_layer_call_and_return_conditional_losses_5670481

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:2*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€V
SigmoidSigmoidBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€Z
IdentityIdentitySigmoid:y:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€2
 
_user_specified_nameinputs
Ѕ
t
H__inference_concatenate_layer_call_and_return_conditional_losses_5671389
inputs_0
inputs_1
identityM
concat/axisConst*
_output_shapes
: *
dtype0*
value	B :w
concatConcatV2inputs_0inputs_1concat/axis:output:0*
N*
T0*'
_output_shapes
:€€€€€€€€€W
IdentityIdentityconcat:output:0*
T0*'
_output_shapes
:€€€€€€€€€"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*9
_input_shapes(
&:€€€€€€€€€:€€€€€€€€€:Q M
'
_output_shapes
:€€€€€€€€€
"
_user_specified_name
inputs_0:QM
'
_output_shapes
:€€€€€€€€€
"
_user_specified_name
inputs_1
≠
L
$__inference__update_step_xla_5671292
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
“
b
I__inference_grad_reverse_layer_call_and_return_conditional_losses_5670374
x

identity_2I
IdentityIdentityx*
T0*'
_output_shapes
:€€€€€€€€€[

Identity_1IdentityIdentity:output:0*
T0*'
_output_shapes
:€€€€€€€€€§
	IdentityN	IdentityNIdentity:output:0x*
T
2*-
_gradient_op_typeCustomGradient-5670368*:
_output_shapes(
&:€€€€€€€€€:€€€€€€€€€\

Identity_2IdentityIdentityN:output:0*
T0*'
_output_shapes
:€€€€€€€€€"!

identity_2Identity_2:output:0*(
_construction_contextkEagerRuntime*&
_input_shapes
:€€€€€€€€€:J F
'
_output_shapes
:€€€€€€€€€

_user_specified_namex
Ы

х
D__inference_dense_5_layer_call_and_return_conditional_losses_5670396

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:2
identityИҐBiasAdd/ReadVariableOpҐMatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:2*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:€€€€€€€€€2P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:€€€€€€€€€2a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:€€€€€€€€€2w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:€€€€€€€€€: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:€€€€€€€€€
 
_user_specified_nameinputs>
$__inference_internal_grad_fn_5671647CustomGradient-5671370>
$__inference_internal_grad_fn_5671656CustomGradient-5670368>
$__inference_internal_grad_fn_5671665CustomGradient-5671161>
$__inference_internal_grad_fn_5671674CustomGradient-5671085>
$__inference_internal_grad_fn_5671683CustomGradient-5670254"Ж
L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*Ь
serving_defaultИ
;
input_10
serving_default_input_1:0€€€€€€€€€	
;
input_20
serving_default_input_2:0€€€€€€€€€7
Adv0
StatefulPartitionedCall:0€€€€€€€€€7
Clf0
StatefulPartitionedCall:1€€€€€€€€€tensorflow/serving/predict:ўЎ
Г
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer-4
layer-5
layer-6
layer_with_weights-3
layer-7
	layer_with_weights-4
	layer-8

layer_with_weights-5

layer-9
layer_with_weights-6
layer-10
layer_with_weights-7
layer-11
layer_with_weights-8
layer-12
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses
_default_save_signature
	optimizer
loss

signatures"
_tf_keras_network
"
_tf_keras_input_layer
ї
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses

kernel
bias"
_tf_keras_layer
ї
 	variables
!trainable_variables
"regularization_losses
#	keras_api
$__call__
*%&call_and_return_all_conditional_losses

&kernel
'bias"
_tf_keras_layer
ї
(	variables
)trainable_variables
*regularization_losses
+	keras_api
,__call__
*-&call_and_return_all_conditional_losses

.kernel
/bias"
_tf_keras_layer
•
0	variables
1trainable_variables
2regularization_losses
3	keras_api
4__call__
*5&call_and_return_all_conditional_losses"
_tf_keras_layer
"
_tf_keras_input_layer
•
6	variables
7trainable_variables
8regularization_losses
9	keras_api
:__call__
*;&call_and_return_all_conditional_losses"
_tf_keras_layer
ї
<	variables
=trainable_variables
>regularization_losses
?	keras_api
@__call__
*A&call_and_return_all_conditional_losses

Bkernel
Cbias"
_tf_keras_layer
ї
D	variables
Etrainable_variables
Fregularization_losses
G	keras_api
H__call__
*I&call_and_return_all_conditional_losses

Jkernel
Kbias"
_tf_keras_layer
ї
L	variables
Mtrainable_variables
Nregularization_losses
O	keras_api
P__call__
*Q&call_and_return_all_conditional_losses

Rkernel
Sbias"
_tf_keras_layer
ї
T	variables
Utrainable_variables
Vregularization_losses
W	keras_api
X__call__
*Y&call_and_return_all_conditional_losses

Zkernel
[bias"
_tf_keras_layer
ї
\	variables
]trainable_variables
^regularization_losses
_	keras_api
`__call__
*a&call_and_return_all_conditional_losses

bkernel
cbias"
_tf_keras_layer
ї
d	variables
etrainable_variables
fregularization_losses
g	keras_api
h__call__
*i&call_and_return_all_conditional_losses

jkernel
kbias"
_tf_keras_layer
¶
0
1
&2
'3
.4
/5
B6
C7
J8
K9
R10
S11
Z12
[13
b14
c15
j16
k17"
trackable_list_wrapper
¶
0
1
&2
'3
.4
/5
B6
C7
J8
K9
R10
S11
Z12
[13
b14
c15
j16
k17"
trackable_list_wrapper
 "
trackable_list_wrapper
 
lnon_trainable_variables

mlayers
nmetrics
olayer_regularization_losses
player_metrics
	variables
trainable_variables
regularization_losses
__call__
_default_save_signature
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses"
_generic_user_object
ў
qtrace_0
rtrace_1
strace_2
ttrace_32о
)__inference_model_1_layer_call_fn_5670530
)__inference_model_1_layer_call_fn_5671016
)__inference_model_1_layer_call_fn_5671060
)__inference_model_1_layer_call_fn_5670818њ
ґ≤≤
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
annotations™ *
 zqtrace_0zrtrace_1zstrace_2zttrace_3
≈
utrace_0
vtrace_1
wtrace_2
xtrace_32Џ
D__inference_model_1_layer_call_and_return_conditional_losses_5671136
D__inference_model_1_layer_call_and_return_conditional_losses_5671212
D__inference_model_1_layer_call_and_return_conditional_losses_5670871
D__inference_model_1_layer_call_and_return_conditional_losses_5670924њ
ґ≤≤
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
annotations™ *
 zutrace_0zvtrace_1zwtrace_2zxtrace_3
÷B”
"__inference__wrapped_model_5670305input_1input_2"Ш
С≤Н
FullArgSpec
argsЪ 
varargsjargs
varkwjkwargs
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotations™ *
 
Ь
y
_variables
z_iterations
{_learning_rate
|_index_dict
}
_momentums
~_velocities
_update_step_xla"
experimentalOptimizer
 "
trackable_list_wrapper
-
Аserving_default"
signature_map
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
≤
Бnon_trainable_variables
Вlayers
Гmetrics
 Дlayer_regularization_losses
Еlayer_metrics
	variables
trainable_variables
regularization_losses
__call__
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses"
_generic_user_object
п
Жtrace_02–
)__inference_dense_3_layer_call_fn_5671311Ґ
Щ≤Х
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
annotations™ *
 zЖtrace_0
К
Зtrace_02л
D__inference_dense_3_layer_call_and_return_conditional_losses_5671322Ґ
Щ≤Х
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
annotations™ *
 zЗtrace_0
 :	22dense_3/kernel
:22dense_3/bias
.
&0
'1"
trackable_list_wrapper
.
&0
'1"
trackable_list_wrapper
 "
trackable_list_wrapper
≤
Иnon_trainable_variables
Йlayers
Кmetrics
 Лlayer_regularization_losses
Мlayer_metrics
 	variables
!trainable_variables
"regularization_losses
$__call__
*%&call_and_return_all_conditional_losses
&%"call_and_return_conditional_losses"
_generic_user_object
п
Нtrace_02–
)__inference_dense_4_layer_call_fn_5671331Ґ
Щ≤Х
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
annotations™ *
 zНtrace_0
К
Оtrace_02л
D__inference_dense_4_layer_call_and_return_conditional_losses_5671342Ґ
Щ≤Х
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
annotations™ *
 zОtrace_0
 :222dense_4/kernel
:22dense_4/bias
.
.0
/1"
trackable_list_wrapper
.
.0
/1"
trackable_list_wrapper
 "
trackable_list_wrapper
≤
Пnon_trainable_variables
Рlayers
Сmetrics
 Тlayer_regularization_losses
Уlayer_metrics
(	variables
)trainable_variables
*regularization_losses
,__call__
*-&call_and_return_all_conditional_losses
&-"call_and_return_conditional_losses"
_generic_user_object
л
Фtrace_02ћ
%__inference_Clf_layer_call_fn_5671351Ґ
Щ≤Х
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
annotations™ *
 zФtrace_0
Ж
Хtrace_02з
@__inference_Clf_layer_call_and_return_conditional_losses_5671362Ґ
Щ≤Х
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
annotations™ *
 zХtrace_0
:22
Clf/kernel
:2Clf/bias
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
≤
Цnon_trainable_variables
Чlayers
Шmetrics
 Щlayer_regularization_losses
Ъlayer_metrics
0	variables
1trainable_variables
2regularization_losses
4__call__
*5&call_and_return_all_conditional_losses
&5"call_and_return_conditional_losses"
_generic_user_object
п
Ыtrace_02–
.__inference_grad_reverse_layer_call_fn_5671367Э
Ф≤Р
FullArgSpec
argsЪ
jself
jx
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotations™ *
 zЫtrace_0
К
Ьtrace_02л
I__inference_grad_reverse_layer_call_and_return_conditional_losses_5671376Э
Ф≤Р
FullArgSpec
argsЪ
jself
jx
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotations™ *
 zЬtrace_0
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
≤
Эnon_trainable_variables
Юlayers
Яmetrics
 †layer_regularization_losses
°layer_metrics
6	variables
7trainable_variables
8regularization_losses
:__call__
*;&call_and_return_all_conditional_losses
&;"call_and_return_conditional_losses"
_generic_user_object
у
Ґtrace_02‘
-__inference_concatenate_layer_call_fn_5671382Ґ
Щ≤Х
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
annotations™ *
 zҐtrace_0
О
£trace_02п
H__inference_concatenate_layer_call_and_return_conditional_losses_5671389Ґ
Щ≤Х
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
annotations™ *
 z£trace_0
.
B0
C1"
trackable_list_wrapper
.
B0
C1"
trackable_list_wrapper
 "
trackable_list_wrapper
≤
§non_trainable_variables
•layers
¶metrics
 Іlayer_regularization_losses
®layer_metrics
<	variables
=trainable_variables
>regularization_losses
@__call__
*A&call_and_return_all_conditional_losses
&A"call_and_return_conditional_losses"
_generic_user_object
п
©trace_02–
)__inference_dense_5_layer_call_fn_5671398Ґ
Щ≤Х
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
annotations™ *
 z©trace_0
К
™trace_02л
D__inference_dense_5_layer_call_and_return_conditional_losses_5671409Ґ
Щ≤Х
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
annotations™ *
 z™trace_0
 :22dense_5/kernel
:22dense_5/bias
.
J0
K1"
trackable_list_wrapper
.
J0
K1"
trackable_list_wrapper
 "
trackable_list_wrapper
≤
Ђnon_trainable_variables
ђlayers
≠metrics
 Ѓlayer_regularization_losses
ѓlayer_metrics
D	variables
Etrainable_variables
Fregularization_losses
H__call__
*I&call_and_return_all_conditional_losses
&I"call_and_return_conditional_losses"
_generic_user_object
п
∞trace_02–
)__inference_dense_6_layer_call_fn_5671418Ґ
Щ≤Х
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
annotations™ *
 z∞trace_0
К
±trace_02л
D__inference_dense_6_layer_call_and_return_conditional_losses_5671429Ґ
Щ≤Х
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
annotations™ *
 z±trace_0
 :222dense_6/kernel
:22dense_6/bias
.
R0
S1"
trackable_list_wrapper
.
R0
S1"
trackable_list_wrapper
 "
trackable_list_wrapper
≤
≤non_trainable_variables
≥layers
іmetrics
 µlayer_regularization_losses
ґlayer_metrics
L	variables
Mtrainable_variables
Nregularization_losses
P__call__
*Q&call_and_return_all_conditional_losses
&Q"call_and_return_conditional_losses"
_generic_user_object
п
Јtrace_02–
)__inference_dense_7_layer_call_fn_5671438Ґ
Щ≤Х
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
annotations™ *
 zЈtrace_0
К
Єtrace_02л
D__inference_dense_7_layer_call_and_return_conditional_losses_5671449Ґ
Щ≤Х
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
annotations™ *
 zЄtrace_0
 :222dense_7/kernel
:22dense_7/bias
.
Z0
[1"
trackable_list_wrapper
.
Z0
[1"
trackable_list_wrapper
 "
trackable_list_wrapper
≤
єnon_trainable_variables
Їlayers
їmetrics
 Љlayer_regularization_losses
љlayer_metrics
T	variables
Utrainable_variables
Vregularization_losses
X__call__
*Y&call_and_return_all_conditional_losses
&Y"call_and_return_conditional_losses"
_generic_user_object
п
Њtrace_02–
)__inference_dense_8_layer_call_fn_5671458Ґ
Щ≤Х
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
annotations™ *
 zЊtrace_0
К
њtrace_02л
D__inference_dense_8_layer_call_and_return_conditional_losses_5671469Ґ
Щ≤Х
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
annotations™ *
 zњtrace_0
 :222dense_8/kernel
:22dense_8/bias
.
b0
c1"
trackable_list_wrapper
.
b0
c1"
trackable_list_wrapper
 "
trackable_list_wrapper
≤
јnon_trainable_variables
Ѕlayers
¬metrics
 √layer_regularization_losses
ƒlayer_metrics
\	variables
]trainable_variables
^regularization_losses
`__call__
*a&call_and_return_all_conditional_losses
&a"call_and_return_conditional_losses"
_generic_user_object
п
≈trace_02–
)__inference_dense_9_layer_call_fn_5671478Ґ
Щ≤Х
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
annotations™ *
 z≈trace_0
К
∆trace_02л
D__inference_dense_9_layer_call_and_return_conditional_losses_5671489Ґ
Щ≤Х
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
annotations™ *
 z∆trace_0
 :222dense_9/kernel
:22dense_9/bias
.
j0
k1"
trackable_list_wrapper
.
j0
k1"
trackable_list_wrapper
 "
trackable_list_wrapper
≤
«non_trainable_variables
»layers
…metrics
  layer_regularization_losses
Ћlayer_metrics
d	variables
etrainable_variables
fregularization_losses
h__call__
*i&call_and_return_all_conditional_losses
&i"call_and_return_conditional_losses"
_generic_user_object
л
ћtrace_02ћ
%__inference_Adv_layer_call_fn_5671498Ґ
Щ≤Х
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
annotations™ *
 zћtrace_0
Ж
Ќtrace_02з
@__inference_Adv_layer_call_and_return_conditional_losses_5671509Ґ
Щ≤Х
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
annotations™ *
 zЌtrace_0
:22
Adv/kernel
:2Adv/bias
 "
trackable_list_wrapper
~
0
1
2
3
4
5
6
7
	8

9
10
11
12"
trackable_list_wrapper
8
ќ0
ѕ1
–2"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
ДBБ
)__inference_model_1_layer_call_fn_5670530input_1input_2"њ
ґ≤≤
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
annotations™ *
 
ЖBГ
)__inference_model_1_layer_call_fn_5671016inputs_0inputs_1"њ
ґ≤≤
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
annotations™ *
 
ЖBГ
)__inference_model_1_layer_call_fn_5671060inputs_0inputs_1"њ
ґ≤≤
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
annotations™ *
 
ДBБ
)__inference_model_1_layer_call_fn_5670818input_1input_2"њ
ґ≤≤
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
annotations™ *
 
°BЮ
D__inference_model_1_layer_call_and_return_conditional_losses_5671136inputs_0inputs_1"њ
ґ≤≤
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
annotations™ *
 
°BЮ
D__inference_model_1_layer_call_and_return_conditional_losses_5671212inputs_0inputs_1"њ
ґ≤≤
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
annotations™ *
 
ЯBЬ
D__inference_model_1_layer_call_and_return_conditional_losses_5670871input_1input_2"њ
ґ≤≤
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
annotations™ *
 
ЯBЬ
D__inference_model_1_layer_call_and_return_conditional_losses_5670924input_1input_2"њ
ґ≤≤
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
annotations™ *
 
в
z0
—1
“2
”3
‘4
’5
÷6
„7
Ў8
ў9
Џ10
џ11
№12
Ё13
ё14
я15
а16
б17
в18
г19
д20
е21
ж22
з23
и24
й25
к26
л27
м28
н29
о30
п31
р32
с33
т34
у35
ф36"
trackable_list_wrapper
:	 2	iteration
: 2learning_rate
 "
trackable_dict_wrapper
Є
—0
”1
’2
„3
ў4
џ5
Ё6
я7
б8
г9
е10
з11
й12
л13
н14
п15
с16
у17"
trackable_list_wrapper
Є
“0
‘1
÷2
Ў3
Џ4
№5
ё6
а7
в8
д9
ж10
и11
к12
м13
о14
р15
т16
ф17"
trackable_list_wrapper
у

хtrace_0
цtrace_1
чtrace_2
шtrace_3
щtrace_4
ъtrace_5
ыtrace_6
ьtrace_7
эtrace_8
юtrace_9
€trace_10
Аtrace_11
Бtrace_12
Вtrace_13
Гtrace_14
Дtrace_15
Еtrace_16
Жtrace_172и
$__inference__update_step_xla_5671217
$__inference__update_step_xla_5671222
$__inference__update_step_xla_5671227
$__inference__update_step_xla_5671232
$__inference__update_step_xla_5671237
$__inference__update_step_xla_5671242
$__inference__update_step_xla_5671247
$__inference__update_step_xla_5671252
$__inference__update_step_xla_5671257
$__inference__update_step_xla_5671262
$__inference__update_step_xla_5671267
$__inference__update_step_xla_5671272
$__inference__update_step_xla_5671277
$__inference__update_step_xla_5671282
$__inference__update_step_xla_5671287
$__inference__update_step_xla_5671292
$__inference__update_step_xla_5671297
$__inference__update_step_xla_5671302є
Ѓ≤™
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
annotations™ *
 0zхtrace_0zцtrace_1zчtrace_2zшtrace_3zщtrace_4zъtrace_5zыtrace_6zьtrace_7zэtrace_8zюtrace_9z€trace_10zАtrace_11zБtrace_12zВtrace_13zГtrace_14zДtrace_15zЕtrace_16zЖtrace_17
”B–
%__inference_signature_wrapper_5670972input_1input_2"Ф
Н≤Й
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
annotations™ *
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
ЁBЏ
)__inference_dense_3_layer_call_fn_5671311inputs"Ґ
Щ≤Х
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
annotations™ *
 
шBх
D__inference_dense_3_layer_call_and_return_conditional_losses_5671322inputs"Ґ
Щ≤Х
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
annotations™ *
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
ЁBЏ
)__inference_dense_4_layer_call_fn_5671331inputs"Ґ
Щ≤Х
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
annotations™ *
 
шBх
D__inference_dense_4_layer_call_and_return_conditional_losses_5671342inputs"Ґ
Щ≤Х
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
annotations™ *
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
ўB÷
%__inference_Clf_layer_call_fn_5671351inputs"Ґ
Щ≤Х
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
annotations™ *
 
фBс
@__inference_Clf_layer_call_and_return_conditional_losses_5671362inputs"Ґ
Щ≤Х
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
annotations™ *
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
ЎB’
.__inference_grad_reverse_layer_call_fn_5671367x"Э
Ф≤Р
FullArgSpec
argsЪ
jself
jx
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotations™ *
 
уBр
I__inference_grad_reverse_layer_call_and_return_conditional_losses_5671376x"Э
Ф≤Р
FullArgSpec
argsЪ
jself
jx
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotations™ *
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
-__inference_concatenate_layer_call_fn_5671382inputs_0inputs_1"Ґ
Щ≤Х
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
annotations™ *
 
ИBЕ
H__inference_concatenate_layer_call_and_return_conditional_losses_5671389inputs_0inputs_1"Ґ
Щ≤Х
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
annotations™ *
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
ЁBЏ
)__inference_dense_5_layer_call_fn_5671398inputs"Ґ
Щ≤Х
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
annotations™ *
 
шBх
D__inference_dense_5_layer_call_and_return_conditional_losses_5671409inputs"Ґ
Щ≤Х
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
annotations™ *
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
ЁBЏ
)__inference_dense_6_layer_call_fn_5671418inputs"Ґ
Щ≤Х
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
annotations™ *
 
шBх
D__inference_dense_6_layer_call_and_return_conditional_losses_5671429inputs"Ґ
Щ≤Х
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
annotations™ *
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
ЁBЏ
)__inference_dense_7_layer_call_fn_5671438inputs"Ґ
Щ≤Х
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
annotations™ *
 
шBх
D__inference_dense_7_layer_call_and_return_conditional_losses_5671449inputs"Ґ
Щ≤Х
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
annotations™ *
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
ЁBЏ
)__inference_dense_8_layer_call_fn_5671458inputs"Ґ
Щ≤Х
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
annotations™ *
 
шBх
D__inference_dense_8_layer_call_and_return_conditional_losses_5671469inputs"Ґ
Щ≤Х
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
annotations™ *
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
ЁBЏ
)__inference_dense_9_layer_call_fn_5671478inputs"Ґ
Щ≤Х
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
annotations™ *
 
шBх
D__inference_dense_9_layer_call_and_return_conditional_losses_5671489inputs"Ґ
Щ≤Х
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
annotations™ *
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
ўB÷
%__inference_Adv_layer_call_fn_5671498inputs"Ґ
Щ≤Х
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
annotations™ *
 
фBс
@__inference_Adv_layer_call_and_return_conditional_losses_5671509inputs"Ґ
Щ≤Х
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
annotations™ *
 
R
З	variables
И	keras_api

Йtotal

Кcount"
_tf_keras_metric
R
Л	variables
М	keras_api

Нtotal

Оcount"
_tf_keras_metric
R
П	variables
Р	keras_api

Сtotal

Тcount"
_tf_keras_metric
%:#	22Adam/m/dense_3/kernel
%:#	22Adam/v/dense_3/kernel
:22Adam/m/dense_3/bias
:22Adam/v/dense_3/bias
%:#222Adam/m/dense_4/kernel
%:#222Adam/v/dense_4/kernel
:22Adam/m/dense_4/bias
:22Adam/v/dense_4/bias
!:22Adam/m/Clf/kernel
!:22Adam/v/Clf/kernel
:2Adam/m/Clf/bias
:2Adam/v/Clf/bias
%:#22Adam/m/dense_5/kernel
%:#22Adam/v/dense_5/kernel
:22Adam/m/dense_5/bias
:22Adam/v/dense_5/bias
%:#222Adam/m/dense_6/kernel
%:#222Adam/v/dense_6/kernel
:22Adam/m/dense_6/bias
:22Adam/v/dense_6/bias
%:#222Adam/m/dense_7/kernel
%:#222Adam/v/dense_7/kernel
:22Adam/m/dense_7/bias
:22Adam/v/dense_7/bias
%:#222Adam/m/dense_8/kernel
%:#222Adam/v/dense_8/kernel
:22Adam/m/dense_8/bias
:22Adam/v/dense_8/bias
%:#222Adam/m/dense_9/kernel
%:#222Adam/v/dense_9/kernel
:22Adam/m/dense_9/bias
:22Adam/v/dense_9/bias
!:22Adam/m/Adv/kernel
!:22Adam/v/Adv/kernel
:2Adam/m/Adv/bias
:2Adam/v/Adv/bias
щBц
$__inference__update_step_xla_5671217gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671222gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671227gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671232gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671237gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671242gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671247gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671252gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671257gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671262gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671267gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671272gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671277gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671282gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671287gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671292gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671297gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
щBц
$__inference__update_step_xla_5671302gradientvariable"Ј
Ѓ≤™
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
annotations™ *
 
0
Й0
К1"
trackable_list_wrapper
.
З	variables"
_generic_user_object
:  (2total
:  (2count
0
Н0
О1"
trackable_list_wrapper
.
Л	variables"
_generic_user_object
:  (2total
:  (2count
0
С0
Т1"
trackable_list_wrapper
.
П	variables"
_generic_user_object
:  (2total
:  (2countІ
@__inference_Adv_layer_call_and_return_conditional_losses_5671509cjk/Ґ,
%Ґ"
 К
inputs€€€€€€€€€2
™ ",Ґ)
"К
tensor_0€€€€€€€€€
Ъ Б
%__inference_Adv_layer_call_fn_5671498Xjk/Ґ,
%Ґ"
 К
inputs€€€€€€€€€2
™ "!К
unknown€€€€€€€€€І
@__inference_Clf_layer_call_and_return_conditional_losses_5671362c.//Ґ,
%Ґ"
 К
inputs€€€€€€€€€2
™ ",Ґ)
"К
tensor_0€€€€€€€€€
Ъ Б
%__inference_Clf_layer_call_fn_5671351X.//Ґ,
%Ґ"
 К
inputs€€€€€€€€€2
™ "!К
unknown€€€€€€€€€Ц
$__inference__update_step_xla_5671217nhҐe
^Ґ[
К
gradient	2
4Т1	Ґ
ъ	2
А
p
` VariableSpec 
`АчєПэ∆?
™ "
 О
$__inference__update_step_xla_5671222f`Ґ]
VҐS
К
gradient2
0Т-	Ґ
ъ2
А
p
` VariableSpec 
`АлоБЕ«?
™ "
 Ц
$__inference__update_step_xla_5671227nhҐe
^Ґ[
К
gradient22
4Т1	Ґ
ъ22
А
p
` VariableSpec 
`јьЊОэ∆?
™ "
 О
$__inference__update_step_xla_5671232f`Ґ]
VҐS
К
gradient2
0Т-	Ґ
ъ2
А
p
` VariableSpec 
`јДЏББ«?
™ "
 Ц
$__inference__update_step_xla_5671237nhҐe
^Ґ[
К
gradient2
4Т1	Ґ
ъ2
А
p
` VariableSpec 
`АОџББ«?
™ "
 О
$__inference__update_step_xla_5671242f`Ґ]
VҐS
К
gradient
0Т-	Ґ
ъ
А
p
` VariableSpec 
`ј‘§Оэ∆?
™ "
 Ц
$__inference__update_step_xla_5671247nhҐe
^Ґ[
К
gradient2
4Т1	Ґ
ъ2
А
p
` VariableSpec 
`†¶ИБЕ«?
™ "
 О
$__inference__update_step_xla_5671252f`Ґ]
VҐS
К
gradient2
0Т-	Ґ
ъ2
А
p
` VariableSpec 
`јЙ§Оэ∆?
™ "
 Ц
$__inference__update_step_xla_5671257nhҐe
^Ґ[
К
gradient22
4Т1	Ґ
ъ22
А
p
` VariableSpec 
`†µђПэ∆?
™ "
 О
$__inference__update_step_xla_5671262f`Ґ]
VҐS
К
gradient2
0Т-	Ґ
ъ2
А
p
` VariableSpec 
`†ИтВБ«?
™ "
 Ц
$__inference__update_step_xla_5671267nhҐe
^Ґ[
К
gradient22
4Т1	Ґ
ъ22
А
p
` VariableSpec 
`†ЉІПэ∆?
™ "
 О
$__inference__update_step_xla_5671272f`Ґ]
VҐS
К
gradient2
0Т-	Ґ
ъ2
А
p
` VariableSpec 
`ј—≠Пэ∆?
™ "
 Ц
$__inference__update_step_xla_5671277nhҐe
^Ґ[
К
gradient22
4Т1	Ґ
ъ22
А
p
` VariableSpec 
`јШ§Оэ∆?
™ "
 О
$__inference__update_step_xla_5671282f`Ґ]
VҐS
К
gradient2
0Т-	Ґ
ъ2
А
p
` VariableSpec 
`†мдВэ∆?
™ "
 Ц
$__inference__update_step_xla_5671287nhҐe
^Ґ[
К
gradient22
4Т1	Ґ
ъ22
А
p
` VariableSpec 
`аЦуВБ«?
™ "
 О
$__inference__update_step_xla_5671292f`Ґ]
VҐS
К
gradient2
0Т-	Ґ
ъ2
А
p
` VariableSpec 
`ј©уВБ«?
™ "
 Ц
$__inference__update_step_xla_5671297nhҐe
^Ґ[
К
gradient2
4Т1	Ґ
ъ2
А
p
` VariableSpec 
`јЯ•Оэ∆?
™ "
 О
$__inference__update_step_xla_5671302f`Ґ]
VҐS
К
gradient
0Т-	Ґ
ъ
А
p
` VariableSpec 
`А“§Оэ∆?
™ "
 ж
"__inference__wrapped_model_5670305њ&'./BCJKRSZ[bcjkXҐU
NҐK
IЪF
!К
input_1€€€€€€€€€	
!К
input_2€€€€€€€€€
™ "O™L
$
AdvК
adv€€€€€€€€€
$
ClfК
clf€€€€€€€€€„
H__inference_concatenate_layer_call_and_return_conditional_losses_5671389КZҐW
PҐM
KЪH
"К
inputs_0€€€€€€€€€
"К
inputs_1€€€€€€€€€
™ ",Ґ)
"К
tensor_0€€€€€€€€€
Ъ ∞
-__inference_concatenate_layer_call_fn_5671382ZҐW
PҐM
KЪH
"К
inputs_0€€€€€€€€€
"К
inputs_1€€€€€€€€€
™ "!К
unknown€€€€€€€€€Ђ
D__inference_dense_3_layer_call_and_return_conditional_losses_5671322c/Ґ,
%Ґ"
 К
inputs€€€€€€€€€	
™ ",Ґ)
"К
tensor_0€€€€€€€€€2
Ъ Е
)__inference_dense_3_layer_call_fn_5671311X/Ґ,
%Ґ"
 К
inputs€€€€€€€€€	
™ "!К
unknown€€€€€€€€€2Ђ
D__inference_dense_4_layer_call_and_return_conditional_losses_5671342c&'/Ґ,
%Ґ"
 К
inputs€€€€€€€€€2
™ ",Ґ)
"К
tensor_0€€€€€€€€€2
Ъ Е
)__inference_dense_4_layer_call_fn_5671331X&'/Ґ,
%Ґ"
 К
inputs€€€€€€€€€2
™ "!К
unknown€€€€€€€€€2Ђ
D__inference_dense_5_layer_call_and_return_conditional_losses_5671409cBC/Ґ,
%Ґ"
 К
inputs€€€€€€€€€
™ ",Ґ)
"К
tensor_0€€€€€€€€€2
Ъ Е
)__inference_dense_5_layer_call_fn_5671398XBC/Ґ,
%Ґ"
 К
inputs€€€€€€€€€
™ "!К
unknown€€€€€€€€€2Ђ
D__inference_dense_6_layer_call_and_return_conditional_losses_5671429cJK/Ґ,
%Ґ"
 К
inputs€€€€€€€€€2
™ ",Ґ)
"К
tensor_0€€€€€€€€€2
Ъ Е
)__inference_dense_6_layer_call_fn_5671418XJK/Ґ,
%Ґ"
 К
inputs€€€€€€€€€2
™ "!К
unknown€€€€€€€€€2Ђ
D__inference_dense_7_layer_call_and_return_conditional_losses_5671449cRS/Ґ,
%Ґ"
 К
inputs€€€€€€€€€2
™ ",Ґ)
"К
tensor_0€€€€€€€€€2
Ъ Е
)__inference_dense_7_layer_call_fn_5671438XRS/Ґ,
%Ґ"
 К
inputs€€€€€€€€€2
™ "!К
unknown€€€€€€€€€2Ђ
D__inference_dense_8_layer_call_and_return_conditional_losses_5671469cZ[/Ґ,
%Ґ"
 К
inputs€€€€€€€€€2
™ ",Ґ)
"К
tensor_0€€€€€€€€€2
Ъ Е
)__inference_dense_8_layer_call_fn_5671458XZ[/Ґ,
%Ґ"
 К
inputs€€€€€€€€€2
™ "!К
unknown€€€€€€€€€2Ђ
D__inference_dense_9_layer_call_and_return_conditional_losses_5671489cbc/Ґ,
%Ґ"
 К
inputs€€€€€€€€€2
™ ",Ґ)
"К
tensor_0€€€€€€€€€2
Ъ Е
)__inference_dense_9_layer_call_fn_5671478Xbc/Ґ,
%Ґ"
 К
inputs€€€€€€€€€2
™ "!К
unknown€€€€€€€€€2І
I__inference_grad_reverse_layer_call_and_return_conditional_losses_5671376Z*Ґ'
 Ґ
К
x€€€€€€€€€
™ ",Ґ)
"К
tensor_0€€€€€€€€€
Ъ Б
.__inference_grad_reverse_layer_call_fn_5671367O*Ґ'
 Ґ
К
x€€€€€€€€€
™ "!К
unknown€€€€€€€€€љ
$__inference_internal_grad_fn_5671647ФeҐb
[ҐX

 
(К%
result_grads_0€€€€€€€€€
(К%
result_grads_1€€€€€€€€€
™ "+Ъ(

 
"К
tensor_1€€€€€€€€€љ
$__inference_internal_grad_fn_5671656ФeҐb
[ҐX

 
(К%
result_grads_0€€€€€€€€€
(К%
result_grads_1€€€€€€€€€
™ "+Ъ(

 
"К
tensor_1€€€€€€€€€љ
$__inference_internal_grad_fn_5671665ФeҐb
[ҐX

 
(К%
result_grads_0€€€€€€€€€
(К%
result_grads_1€€€€€€€€€
™ "+Ъ(

 
"К
tensor_1€€€€€€€€€љ
$__inference_internal_grad_fn_5671674ФeҐb
[ҐX

 
(К%
result_grads_0€€€€€€€€€
(К%
result_grads_1€€€€€€€€€
™ "+Ъ(

 
"К
tensor_1€€€€€€€€€љ
$__inference_internal_grad_fn_5671683ФeҐb
[ҐX

 
(К%
result_grads_0€€€€€€€€€
(К%
result_grads_1€€€€€€€€€
™ "+Ъ(

 
"К
tensor_1€€€€€€€€€Ъ
D__inference_model_1_layer_call_and_return_conditional_losses_5670871—&'./BCJKRSZ[bcjk`Ґ]
VҐS
IЪF
!К
input_1€€€€€€€€€	
!К
input_2€€€€€€€€€
p 

 
™ "YҐV
OЪL
$К!

tensor_0_0€€€€€€€€€
$К!

tensor_0_1€€€€€€€€€
Ъ Ъ
D__inference_model_1_layer_call_and_return_conditional_losses_5670924—&'./BCJKRSZ[bcjk`Ґ]
VҐS
IЪF
!К
input_1€€€€€€€€€	
!К
input_2€€€€€€€€€
p

 
™ "YҐV
OЪL
$К!

tensor_0_0€€€€€€€€€
$К!

tensor_0_1€€€€€€€€€
Ъ Ь
D__inference_model_1_layer_call_and_return_conditional_losses_5671136”&'./BCJKRSZ[bcjkbҐ_
XҐU
KЪH
"К
inputs_0€€€€€€€€€	
"К
inputs_1€€€€€€€€€
p 

 
™ "YҐV
OЪL
$К!

tensor_0_0€€€€€€€€€
$К!

tensor_0_1€€€€€€€€€
Ъ Ь
D__inference_model_1_layer_call_and_return_conditional_losses_5671212”&'./BCJKRSZ[bcjkbҐ_
XҐU
KЪH
"К
inputs_0€€€€€€€€€	
"К
inputs_1€€€€€€€€€
p

 
™ "YҐV
OЪL
$К!

tensor_0_0€€€€€€€€€
$К!

tensor_0_1€€€€€€€€€
Ъ с
)__inference_model_1_layer_call_fn_5670530√&'./BCJKRSZ[bcjk`Ґ]
VҐS
IЪF
!К
input_1€€€€€€€€€	
!К
input_2€€€€€€€€€
p 

 
™ "KЪH
"К
tensor_0€€€€€€€€€
"К
tensor_1€€€€€€€€€с
)__inference_model_1_layer_call_fn_5670818√&'./BCJKRSZ[bcjk`Ґ]
VҐS
IЪF
!К
input_1€€€€€€€€€	
!К
input_2€€€€€€€€€
p

 
™ "KЪH
"К
tensor_0€€€€€€€€€
"К
tensor_1€€€€€€€€€у
)__inference_model_1_layer_call_fn_5671016≈&'./BCJKRSZ[bcjkbҐ_
XҐU
KЪH
"К
inputs_0€€€€€€€€€	
"К
inputs_1€€€€€€€€€
p 

 
™ "KЪH
"К
tensor_0€€€€€€€€€
"К
tensor_1€€€€€€€€€у
)__inference_model_1_layer_call_fn_5671060≈&'./BCJKRSZ[bcjkbҐ_
XҐU
KЪH
"К
inputs_0€€€€€€€€€	
"К
inputs_1€€€€€€€€€
p

 
™ "KЪH
"К
tensor_0€€€€€€€€€
"К
tensor_1€€€€€€€€€ъ
%__inference_signature_wrapper_5670972–&'./BCJKRSZ[bcjkiҐf
Ґ 
_™\
,
input_1!К
input_1€€€€€€€€€	
,
input_2!К
input_2€€€€€€€€€"O™L
$
AdvК
adv€€€€€€€€€
$
ClfК
clf€€€€€€€€€