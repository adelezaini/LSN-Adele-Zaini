��
��
8
Const
output"dtype"
valuetensor"
dtypetype

NoOp
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetype�
�
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
executor_typestring �
q
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape�"serve*2.0.02unknown8��
|
dense_600/kernelVarHandleOp*!
shared_namedense_600/kernel*
dtype0*
_output_shapes
: *
shape
:
u
$dense_600/kernel/Read/ReadVariableOpReadVariableOpdense_600/kernel*
_output_shapes

:*
dtype0
t
dense_600/biasVarHandleOp*
shape:*
shared_namedense_600/bias*
dtype0*
_output_shapes
: 
m
"dense_600/bias/Read/ReadVariableOpReadVariableOpdense_600/bias*
dtype0*
_output_shapes
:
|
dense_601/kernelVarHandleOp*
shape
:*!
shared_namedense_601/kernel*
dtype0*
_output_shapes
: 
u
$dense_601/kernel/Read/ReadVariableOpReadVariableOpdense_601/kernel*
dtype0*
_output_shapes

:
t
dense_601/biasVarHandleOp*
shape:*
shared_namedense_601/bias*
dtype0*
_output_shapes
: 
m
"dense_601/bias/Read/ReadVariableOpReadVariableOpdense_601/bias*
dtype0*
_output_shapes
:
|
dense_602/kernelVarHandleOp*
dtype0*
_output_shapes
: *
shape
:*!
shared_namedense_602/kernel
u
$dense_602/kernel/Read/ReadVariableOpReadVariableOpdense_602/kernel*
dtype0*
_output_shapes

:
t
dense_602/biasVarHandleOp*
dtype0*
_output_shapes
: *
shape:*
shared_namedense_602/bias
m
"dense_602/bias/Read/ReadVariableOpReadVariableOpdense_602/bias*
dtype0*
_output_shapes
:
|
dense_603/kernelVarHandleOp*
dtype0*
_output_shapes
: *
shape
:*!
shared_namedense_603/kernel
u
$dense_603/kernel/Read/ReadVariableOpReadVariableOpdense_603/kernel*
dtype0*
_output_shapes

:
t
dense_603/biasVarHandleOp*
shape:*
shared_namedense_603/bias*
dtype0*
_output_shapes
: 
m
"dense_603/bias/Read/ReadVariableOpReadVariableOpdense_603/bias*
dtype0*
_output_shapes
:
d
SGD/iterVarHandleOp*
dtype0	*
_output_shapes
: *
shape: *
shared_name
SGD/iter
]
SGD/iter/Read/ReadVariableOpReadVariableOpSGD/iter*
dtype0	*
_output_shapes
: 
f
	SGD/decayVarHandleOp*
shape: *
shared_name	SGD/decay*
dtype0*
_output_shapes
: 
_
SGD/decay/Read/ReadVariableOpReadVariableOp	SGD/decay*
dtype0*
_output_shapes
: 
v
SGD/learning_rateVarHandleOp*"
shared_nameSGD/learning_rate*
dtype0*
_output_shapes
: *
shape: 
o
%SGD/learning_rate/Read/ReadVariableOpReadVariableOpSGD/learning_rate*
dtype0*
_output_shapes
: 
l
SGD/momentumVarHandleOp*
shared_nameSGD/momentum*
dtype0*
_output_shapes
: *
shape: 
e
 SGD/momentum/Read/ReadVariableOpReadVariableOpSGD/momentum*
dtype0*
_output_shapes
: 
^
totalVarHandleOp*
shared_nametotal*
dtype0*
_output_shapes
: *
shape: 
W
total/Read/ReadVariableOpReadVariableOptotal*
dtype0*
_output_shapes
: 
^
countVarHandleOp*
shape: *
shared_namecount*
dtype0*
_output_shapes
: 
W
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
dtype0

NoOpNoOp
�
ConstConst"/device:CPU:0*�
value�B� B�
�
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer_with_weights-3
layer-4
	optimizer
trainable_variables
	variables
	regularization_losses

	keras_api

signatures
R
trainable_variables
	variables
regularization_losses
	keras_api
h

kernel
bias
trainable_variables
	variables
regularization_losses
	keras_api
h

kernel
bias
trainable_variables
	variables
regularization_losses
	keras_api
h

kernel
bias
trainable_variables
	variables
 regularization_losses
!	keras_api
h

"kernel
#bias
$trainable_variables
%	variables
&regularization_losses
'	keras_api
6
(iter
	)decay
*learning_rate
+momentum
8
0
1
2
3
4
5
"6
#7
8
0
1
2
3
4
5
"6
#7
 
�

,layers
-layer_regularization_losses
.non_trainable_variables
trainable_variables
/metrics
	variables
	regularization_losses
 
 
 
 
�

0layers
1layer_regularization_losses
2non_trainable_variables
trainable_variables
3metrics
	variables
regularization_losses
\Z
VARIABLE_VALUEdense_600/kernel6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE
XV
VARIABLE_VALUEdense_600/bias4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE

0
1

0
1
 
�

4layers
5layer_regularization_losses
6non_trainable_variables
trainable_variables
7metrics
	variables
regularization_losses
\Z
VARIABLE_VALUEdense_601/kernel6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE
XV
VARIABLE_VALUEdense_601/bias4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE

0
1

0
1
 
�

8layers
9layer_regularization_losses
:non_trainable_variables
trainable_variables
;metrics
	variables
regularization_losses
\Z
VARIABLE_VALUEdense_602/kernel6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE
XV
VARIABLE_VALUEdense_602/bias4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUE

0
1

0
1
 
�

<layers
=layer_regularization_losses
>non_trainable_variables
trainable_variables
?metrics
	variables
 regularization_losses
\Z
VARIABLE_VALUEdense_603/kernel6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUE
XV
VARIABLE_VALUEdense_603/bias4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUE

"0
#1

"0
#1
 
�

@layers
Alayer_regularization_losses
Bnon_trainable_variables
$trainable_variables
Cmetrics
%	variables
&regularization_losses
GE
VARIABLE_VALUESGD/iter)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUE
IG
VARIABLE_VALUE	SGD/decay*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUE
YW
VARIABLE_VALUESGD/learning_rate2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE
OM
VARIABLE_VALUESGD/momentum-optimizer/momentum/.ATTRIBUTES/VARIABLE_VALUE

0
1
2
3
 
 

D0
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
x
	Etotal
	Fcount
G
_fn_kwargs
Htrainable_variables
I	variables
Jregularization_losses
K	keras_api
OM
VARIABLE_VALUEtotal4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE
OM
VARIABLE_VALUEcount4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE
 
 

E0
F1
 
�

Llayers
Mlayer_regularization_losses
Nnon_trainable_variables
Htrainable_variables
Ometrics
I	variables
Jregularization_losses
 
 

E0
F1
 *
dtype0*
_output_shapes
: 
�
serving_default_dense_600_inputPlaceholder*'
_output_shapes
:���������*
shape:���������*
dtype0
�
StatefulPartitionedCallStatefulPartitionedCallserving_default_dense_600_inputdense_600/kerneldense_600/biasdense_601/kerneldense_601/biasdense_602/kerneldense_602/biasdense_603/kerneldense_603/bias*.
_gradient_op_typePartitionedCall-1123883*.
f)R'
%__inference_signature_wrapper_1123627*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2	*'
_output_shapes
:���������
O
saver_filenamePlaceholder*
shape: *
dtype0*
_output_shapes
: 
�
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename$dense_600/kernel/Read/ReadVariableOp"dense_600/bias/Read/ReadVariableOp$dense_601/kernel/Read/ReadVariableOp"dense_601/bias/Read/ReadVariableOp$dense_602/kernel/Read/ReadVariableOp"dense_602/bias/Read/ReadVariableOp$dense_603/kernel/Read/ReadVariableOp"dense_603/bias/Read/ReadVariableOpSGD/iter/Read/ReadVariableOpSGD/decay/Read/ReadVariableOp%SGD/learning_rate/Read/ReadVariableOp SGD/momentum/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOpConst*.
_gradient_op_typePartitionedCall-1123919*)
f$R"
 __inference__traced_save_1123918*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2	*
_output_shapes
: 
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenamedense_600/kerneldense_600/biasdense_601/kerneldense_601/biasdense_602/kerneldense_602/biasdense_603/kerneldense_603/biasSGD/iter	SGD/decaySGD/learning_rateSGD/momentumtotalcount*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*
_output_shapes
: *.
_gradient_op_typePartitionedCall-1123974*,
f'R%
#__inference__traced_restore_1123973��
�
�
F__inference_dense_600_layer_call_and_return_conditional_losses_1123413

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
_output_shapes

:*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*'
_output_shapes
:���������*
T0�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������N
EluEluBiasAdd:output:0*'
_output_shapes
:���������*
T0N
	Greater/yConst*
_output_shapes
: *
valueB
 *    *
dtype0j
GreaterGreaterBiasAdd:output:0Greater/y:output:0*'
_output_shapes
:���������*
T0J
mul/xConst*
valueB
 *}-�?*
dtype0*
_output_shapes
: _
mulMulmul/x:output:0Elu:activations:0*
T0*'
_output_shapes
:���������k
SelectSelectGreater:z:0Elu:activations:0mul:z:0*
T0*'
_output_shapes
:���������L
mul_1/xConst*
valueB
 *_}�?*
dtype0*
_output_shapes
: a
mul_1Mulmul_1/x:output:0Select:output:0*'
_output_shapes
:���������*
T0�
IdentityIdentity	mul_1:z:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:& "
 
_user_specified_nameinputs: : 
�8
�
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123733

inputs,
(dense_600_matmul_readvariableop_resource-
)dense_600_biasadd_readvariableop_resource,
(dense_601_matmul_readvariableop_resource-
)dense_601_biasadd_readvariableop_resource,
(dense_602_matmul_readvariableop_resource-
)dense_602_biasadd_readvariableop_resource,
(dense_603_matmul_readvariableop_resource-
)dense_603_biasadd_readvariableop_resource
identity�� dense_600/BiasAdd/ReadVariableOp�dense_600/MatMul/ReadVariableOp� dense_601/BiasAdd/ReadVariableOp�dense_601/MatMul/ReadVariableOp� dense_602/BiasAdd/ReadVariableOp�dense_602/MatMul/ReadVariableOp� dense_603/BiasAdd/ReadVariableOp�dense_603/MatMul/ReadVariableOp�
dense_600/MatMul/ReadVariableOpReadVariableOp(dense_600_matmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
_output_shapes

:*
dtype0}
dense_600/MatMulMatMulinputs'dense_600/MatMul/ReadVariableOp:value:0*'
_output_shapes
:���������*
T0�
 dense_600/BiasAdd/ReadVariableOpReadVariableOp)dense_600_biasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:�
dense_600/BiasAddBiasAdddense_600/MatMul:product:0(dense_600/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������b
dense_600/EluEludense_600/BiasAdd:output:0*'
_output_shapes
:���������*
T0X
dense_600/Greater/yConst*
valueB
 *    *
dtype0*
_output_shapes
: �
dense_600/GreaterGreaterdense_600/BiasAdd:output:0dense_600/Greater/y:output:0*
T0*'
_output_shapes
:���������T
dense_600/mul/xConst*
valueB
 *}-�?*
dtype0*
_output_shapes
: }
dense_600/mulMuldense_600/mul/x:output:0dense_600/Elu:activations:0*
T0*'
_output_shapes
:����������
dense_600/SelectSelectdense_600/Greater:z:0dense_600/Elu:activations:0dense_600/mul:z:0*
T0*'
_output_shapes
:���������V
dense_600/mul_1/xConst*
valueB
 *_}�?*
dtype0*
_output_shapes
: 
dense_600/mul_1Muldense_600/mul_1/x:output:0dense_600/Select:output:0*
T0*'
_output_shapes
:����������
dense_601/MatMul/ReadVariableOpReadVariableOp(dense_601_matmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:�
dense_601/MatMulMatMuldense_600/mul_1:z:0'dense_601/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
 dense_601/BiasAdd/ReadVariableOpReadVariableOp)dense_601_biasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:�
dense_601/BiasAddBiasAdddense_601/MatMul:product:0(dense_601/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������b
dense_601/EluEludense_601/BiasAdd:output:0*
T0*'
_output_shapes
:���������X
dense_601/Greater/yConst*
dtype0*
_output_shapes
: *
valueB
 *    �
dense_601/GreaterGreaterdense_601/BiasAdd:output:0dense_601/Greater/y:output:0*
T0*'
_output_shapes
:���������T
dense_601/mul/xConst*
valueB
 *}-�?*
dtype0*
_output_shapes
: }
dense_601/mulMuldense_601/mul/x:output:0dense_601/Elu:activations:0*'
_output_shapes
:���������*
T0�
dense_601/SelectSelectdense_601/Greater:z:0dense_601/Elu:activations:0dense_601/mul:z:0*'
_output_shapes
:���������*
T0V
dense_601/mul_1/xConst*
valueB
 *_}�?*
dtype0*
_output_shapes
: 
dense_601/mul_1Muldense_601/mul_1/x:output:0dense_601/Select:output:0*'
_output_shapes
:���������*
T0�
dense_602/MatMul/ReadVariableOpReadVariableOp(dense_602_matmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:�
dense_602/MatMulMatMuldense_601/mul_1:z:0'dense_602/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
 dense_602/BiasAdd/ReadVariableOpReadVariableOp)dense_602_biasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:�
dense_602/BiasAddBiasAdddense_602/MatMul:product:0(dense_602/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������b
dense_602/EluEludense_602/BiasAdd:output:0*
T0*'
_output_shapes
:���������X
dense_602/Greater/yConst*
valueB
 *    *
dtype0*
_output_shapes
: �
dense_602/GreaterGreaterdense_602/BiasAdd:output:0dense_602/Greater/y:output:0*'
_output_shapes
:���������*
T0T
dense_602/mul/xConst*
valueB
 *}-�?*
dtype0*
_output_shapes
: }
dense_602/mulMuldense_602/mul/x:output:0dense_602/Elu:activations:0*
T0*'
_output_shapes
:����������
dense_602/SelectSelectdense_602/Greater:z:0dense_602/Elu:activations:0dense_602/mul:z:0*'
_output_shapes
:���������*
T0V
dense_602/mul_1/xConst*
valueB
 *_}�?*
dtype0*
_output_shapes
: 
dense_602/mul_1Muldense_602/mul_1/x:output:0dense_602/Select:output:0*
T0*'
_output_shapes
:����������
dense_603/MatMul/ReadVariableOpReadVariableOp(dense_603_matmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:�
dense_603/MatMulMatMuldense_602/mul_1:z:0'dense_603/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
 dense_603/BiasAdd/ReadVariableOpReadVariableOp)dense_603_biasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:�
dense_603/BiasAddBiasAdddense_603/MatMul:product:0(dense_603/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
IdentityIdentitydense_603/BiasAdd:output:0!^dense_600/BiasAdd/ReadVariableOp ^dense_600/MatMul/ReadVariableOp!^dense_601/BiasAdd/ReadVariableOp ^dense_601/MatMul/ReadVariableOp!^dense_602/BiasAdd/ReadVariableOp ^dense_602/MatMul/ReadVariableOp!^dense_603/BiasAdd/ReadVariableOp ^dense_603/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::2B
dense_601/MatMul/ReadVariableOpdense_601/MatMul/ReadVariableOp2B
dense_603/MatMul/ReadVariableOpdense_603/MatMul/ReadVariableOp2D
 dense_603/BiasAdd/ReadVariableOp dense_603/BiasAdd/ReadVariableOp2D
 dense_602/BiasAdd/ReadVariableOp dense_602/BiasAdd/ReadVariableOp2B
dense_600/MatMul/ReadVariableOpdense_600/MatMul/ReadVariableOp2D
 dense_601/BiasAdd/ReadVariableOp dense_601/BiasAdd/ReadVariableOp2B
dense_602/MatMul/ReadVariableOpdense_602/MatMul/ReadVariableOp2D
 dense_600/BiasAdd/ReadVariableOp dense_600/BiasAdd/ReadVariableOp: : : : : : : : :& "
 
_user_specified_nameinputs
�

�
/__inference_sequential_87_layer_call_fn_1123746

inputs"
statefulpartitionedcall_args_1"
statefulpartitionedcall_args_2"
statefulpartitionedcall_args_3"
statefulpartitionedcall_args_4"
statefulpartitionedcall_args_5"
statefulpartitionedcall_args_6"
statefulpartitionedcall_args_7"
statefulpartitionedcall_args_8
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsstatefulpartitionedcall_args_1statefulpartitionedcall_args_2statefulpartitionedcall_args_3statefulpartitionedcall_args_4statefulpartitionedcall_args_5statefulpartitionedcall_args_6statefulpartitionedcall_args_7statefulpartitionedcall_args_8*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2	*'
_output_shapes
:���������*.
_gradient_op_typePartitionedCall-1123566*S
fNRL
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123565�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*'
_output_shapes
:���������*
T0"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::22
StatefulPartitionedCallStatefulPartitionedCall:& "
 
_user_specified_nameinputs: : : : : : : : 
�
�
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123546
dense_600_input,
(dense_600_statefulpartitionedcall_args_1,
(dense_600_statefulpartitionedcall_args_2,
(dense_601_statefulpartitionedcall_args_1,
(dense_601_statefulpartitionedcall_args_2,
(dense_602_statefulpartitionedcall_args_1,
(dense_602_statefulpartitionedcall_args_2,
(dense_603_statefulpartitionedcall_args_1,
(dense_603_statefulpartitionedcall_args_2
identity��!dense_600/StatefulPartitionedCall�!dense_601/StatefulPartitionedCall�!dense_602/StatefulPartitionedCall�!dense_603/StatefulPartitionedCall�
!dense_600/StatefulPartitionedCallStatefulPartitionedCalldense_600_input(dense_600_statefulpartitionedcall_args_1(dense_600_statefulpartitionedcall_args_2*.
_gradient_op_typePartitionedCall-1123419*O
fJRH
F__inference_dense_600_layer_call_and_return_conditional_losses_1123413*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:����������
!dense_601/StatefulPartitionedCallStatefulPartitionedCall*dense_600/StatefulPartitionedCall:output:0(dense_601_statefulpartitionedcall_args_1(dense_601_statefulpartitionedcall_args_2**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:���������*.
_gradient_op_typePartitionedCall-1123454*O
fJRH
F__inference_dense_601_layer_call_and_return_conditional_losses_1123448*
Tout
2�
!dense_602/StatefulPartitionedCallStatefulPartitionedCall*dense_601/StatefulPartitionedCall:output:0(dense_602_statefulpartitionedcall_args_1(dense_602_statefulpartitionedcall_args_2*
Tin
2*'
_output_shapes
:���������*.
_gradient_op_typePartitionedCall-1123489*O
fJRH
F__inference_dense_602_layer_call_and_return_conditional_losses_1123483*
Tout
2**
config_proto

CPU

GPU 2J 8�
!dense_603/StatefulPartitionedCallStatefulPartitionedCall*dense_602/StatefulPartitionedCall:output:0(dense_603_statefulpartitionedcall_args_1(dense_603_statefulpartitionedcall_args_2*O
fJRH
F__inference_dense_603_layer_call_and_return_conditional_losses_1123510*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:���������*.
_gradient_op_typePartitionedCall-1123516�
IdentityIdentity*dense_603/StatefulPartitionedCall:output:0"^dense_600/StatefulPartitionedCall"^dense_601/StatefulPartitionedCall"^dense_602/StatefulPartitionedCall"^dense_603/StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::2F
!dense_600/StatefulPartitionedCall!dense_600/StatefulPartitionedCall2F
!dense_601/StatefulPartitionedCall!dense_601/StatefulPartitionedCall2F
!dense_602/StatefulPartitionedCall!dense_602/StatefulPartitionedCall2F
!dense_603/StatefulPartitionedCall!dense_603/StatefulPartitionedCall:/ +
)
_user_specified_namedense_600_input: : : : : : : : 
�

�
/__inference_sequential_87_layer_call_fn_1123577
dense_600_input"
statefulpartitionedcall_args_1"
statefulpartitionedcall_args_2"
statefulpartitionedcall_args_3"
statefulpartitionedcall_args_4"
statefulpartitionedcall_args_5"
statefulpartitionedcall_args_6"
statefulpartitionedcall_args_7"
statefulpartitionedcall_args_8
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCalldense_600_inputstatefulpartitionedcall_args_1statefulpartitionedcall_args_2statefulpartitionedcall_args_3statefulpartitionedcall_args_4statefulpartitionedcall_args_5statefulpartitionedcall_args_6statefulpartitionedcall_args_7statefulpartitionedcall_args_8*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2	*'
_output_shapes
:���������*.
_gradient_op_typePartitionedCall-1123566*S
fNRL
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123565�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::22
StatefulPartitionedCallStatefulPartitionedCall:/ +
)
_user_specified_namedense_600_input: : : : : : : : 
�

�
/__inference_sequential_87_layer_call_fn_1123609
dense_600_input"
statefulpartitionedcall_args_1"
statefulpartitionedcall_args_2"
statefulpartitionedcall_args_3"
statefulpartitionedcall_args_4"
statefulpartitionedcall_args_5"
statefulpartitionedcall_args_6"
statefulpartitionedcall_args_7"
statefulpartitionedcall_args_8
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCalldense_600_inputstatefulpartitionedcall_args_1statefulpartitionedcall_args_2statefulpartitionedcall_args_3statefulpartitionedcall_args_4statefulpartitionedcall_args_5statefulpartitionedcall_args_6statefulpartitionedcall_args_7statefulpartitionedcall_args_8*'
_output_shapes
:���������*
Tin
2	*.
_gradient_op_typePartitionedCall-1123598*S
fNRL
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123597*
Tout
2**
config_proto

CPU

GPU 2J 8�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::22
StatefulPartitionedCallStatefulPartitionedCall:/ +
)
_user_specified_namedense_600_input: : : : : : : : 
�
�
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123565

inputs,
(dense_600_statefulpartitionedcall_args_1,
(dense_600_statefulpartitionedcall_args_2,
(dense_601_statefulpartitionedcall_args_1,
(dense_601_statefulpartitionedcall_args_2,
(dense_602_statefulpartitionedcall_args_1,
(dense_602_statefulpartitionedcall_args_2,
(dense_603_statefulpartitionedcall_args_1,
(dense_603_statefulpartitionedcall_args_2
identity��!dense_600/StatefulPartitionedCall�!dense_601/StatefulPartitionedCall�!dense_602/StatefulPartitionedCall�!dense_603/StatefulPartitionedCall�
!dense_600/StatefulPartitionedCallStatefulPartitionedCallinputs(dense_600_statefulpartitionedcall_args_1(dense_600_statefulpartitionedcall_args_2*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:���������*.
_gradient_op_typePartitionedCall-1123419*O
fJRH
F__inference_dense_600_layer_call_and_return_conditional_losses_1123413�
!dense_601/StatefulPartitionedCallStatefulPartitionedCall*dense_600/StatefulPartitionedCall:output:0(dense_601_statefulpartitionedcall_args_1(dense_601_statefulpartitionedcall_args_2*
Tin
2*'
_output_shapes
:���������*.
_gradient_op_typePartitionedCall-1123454*O
fJRH
F__inference_dense_601_layer_call_and_return_conditional_losses_1123448*
Tout
2**
config_proto

CPU

GPU 2J 8�
!dense_602/StatefulPartitionedCallStatefulPartitionedCall*dense_601/StatefulPartitionedCall:output:0(dense_602_statefulpartitionedcall_args_1(dense_602_statefulpartitionedcall_args_2**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:���������*.
_gradient_op_typePartitionedCall-1123489*O
fJRH
F__inference_dense_602_layer_call_and_return_conditional_losses_1123483*
Tout
2�
!dense_603/StatefulPartitionedCallStatefulPartitionedCall*dense_602/StatefulPartitionedCall:output:0(dense_603_statefulpartitionedcall_args_1(dense_603_statefulpartitionedcall_args_2*.
_gradient_op_typePartitionedCall-1123516*O
fJRH
F__inference_dense_603_layer_call_and_return_conditional_losses_1123510*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:����������
IdentityIdentity*dense_603/StatefulPartitionedCall:output:0"^dense_600/StatefulPartitionedCall"^dense_601/StatefulPartitionedCall"^dense_602/StatefulPartitionedCall"^dense_603/StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::2F
!dense_600/StatefulPartitionedCall!dense_600/StatefulPartitionedCall2F
!dense_601/StatefulPartitionedCall!dense_601/StatefulPartitionedCall2F
!dense_602/StatefulPartitionedCall!dense_602/StatefulPartitionedCall2F
!dense_603/StatefulPartitionedCall!dense_603/StatefulPartitionedCall: : :& "
 
_user_specified_nameinputs: : : : : : 
�
�
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123528
dense_600_input,
(dense_600_statefulpartitionedcall_args_1,
(dense_600_statefulpartitionedcall_args_2,
(dense_601_statefulpartitionedcall_args_1,
(dense_601_statefulpartitionedcall_args_2,
(dense_602_statefulpartitionedcall_args_1,
(dense_602_statefulpartitionedcall_args_2,
(dense_603_statefulpartitionedcall_args_1,
(dense_603_statefulpartitionedcall_args_2
identity��!dense_600/StatefulPartitionedCall�!dense_601/StatefulPartitionedCall�!dense_602/StatefulPartitionedCall�!dense_603/StatefulPartitionedCall�
!dense_600/StatefulPartitionedCallStatefulPartitionedCalldense_600_input(dense_600_statefulpartitionedcall_args_1(dense_600_statefulpartitionedcall_args_2**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:���������*.
_gradient_op_typePartitionedCall-1123419*O
fJRH
F__inference_dense_600_layer_call_and_return_conditional_losses_1123413*
Tout
2�
!dense_601/StatefulPartitionedCallStatefulPartitionedCall*dense_600/StatefulPartitionedCall:output:0(dense_601_statefulpartitionedcall_args_1(dense_601_statefulpartitionedcall_args_2*.
_gradient_op_typePartitionedCall-1123454*O
fJRH
F__inference_dense_601_layer_call_and_return_conditional_losses_1123448*
Tout
2**
config_proto

CPU

GPU 2J 8*'
_output_shapes
:���������*
Tin
2�
!dense_602/StatefulPartitionedCallStatefulPartitionedCall*dense_601/StatefulPartitionedCall:output:0(dense_602_statefulpartitionedcall_args_1(dense_602_statefulpartitionedcall_args_2**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:���������*.
_gradient_op_typePartitionedCall-1123489*O
fJRH
F__inference_dense_602_layer_call_and_return_conditional_losses_1123483*
Tout
2�
!dense_603/StatefulPartitionedCallStatefulPartitionedCall*dense_602/StatefulPartitionedCall:output:0(dense_603_statefulpartitionedcall_args_1(dense_603_statefulpartitionedcall_args_2**
config_proto

CPU

GPU 2J 8*'
_output_shapes
:���������*
Tin
2*.
_gradient_op_typePartitionedCall-1123516*O
fJRH
F__inference_dense_603_layer_call_and_return_conditional_losses_1123510*
Tout
2�
IdentityIdentity*dense_603/StatefulPartitionedCall:output:0"^dense_600/StatefulPartitionedCall"^dense_601/StatefulPartitionedCall"^dense_602/StatefulPartitionedCall"^dense_603/StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::2F
!dense_603/StatefulPartitionedCall!dense_603/StatefulPartitionedCall2F
!dense_600/StatefulPartitionedCall!dense_600/StatefulPartitionedCall2F
!dense_601/StatefulPartitionedCall!dense_601/StatefulPartitionedCall2F
!dense_602/StatefulPartitionedCall!dense_602/StatefulPartitionedCall: : : : : : :/ +
)
_user_specified_namedense_600_input: : 
�%
�
 __inference__traced_save_1123918
file_prefix/
+savev2_dense_600_kernel_read_readvariableop-
)savev2_dense_600_bias_read_readvariableop/
+savev2_dense_601_kernel_read_readvariableop-
)savev2_dense_601_bias_read_readvariableop/
+savev2_dense_602_kernel_read_readvariableop-
)savev2_dense_602_bias_read_readvariableop/
+savev2_dense_603_kernel_read_readvariableop-
)savev2_dense_603_bias_read_readvariableop'
#savev2_sgd_iter_read_readvariableop	(
$savev2_sgd_decay_read_readvariableop0
,savev2_sgd_learning_rate_read_readvariableop+
'savev2_sgd_momentum_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop
savev2_1_const

identity_1��MergeV2Checkpoints�SaveV2�SaveV2_1�
StringJoin/inputs_1Const"/device:CPU:0*<
value3B1 B+_temp_27ee43e34a7c4a40aee844fa2ae85874/part*
dtype0*
_output_shapes
: s

StringJoin
StringJoinfile_prefixStringJoin/inputs_1:output:0"/device:CPU:0*
_output_shapes
: *
NL

num_shardsConst*
value	B :*
dtype0*
_output_shapes
: f
ShardedFilename/shardConst"/device:CPU:0*
dtype0*
_output_shapes
: *
value	B : �
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: �
SaveV2/tensor_namesConst"/device:CPU:0*�
value�B�B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB-optimizer/momentum/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:�
SaveV2/shape_and_slicesConst"/device:CPU:0*/
value&B$B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:�
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0+savev2_dense_600_kernel_read_readvariableop)savev2_dense_600_bias_read_readvariableop+savev2_dense_601_kernel_read_readvariableop)savev2_dense_601_bias_read_readvariableop+savev2_dense_602_kernel_read_readvariableop)savev2_dense_602_bias_read_readvariableop+savev2_dense_603_kernel_read_readvariableop)savev2_dense_603_bias_read_readvariableop#savev2_sgd_iter_read_readvariableop$savev2_sgd_decay_read_readvariableop,savev2_sgd_learning_rate_read_readvariableop'savev2_sgd_momentum_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableop"/device:CPU:0*
_output_shapes
 *
dtypes
2	h
ShardedFilename_1/shardConst"/device:CPU:0*
value	B :*
dtype0*
_output_shapes
: �
ShardedFilename_1ShardedFilenameStringJoin:output:0 ShardedFilename_1/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: �
SaveV2_1/tensor_namesConst"/device:CPU:0*1
value(B&B_CHECKPOINTABLE_OBJECT_GRAPH*
dtype0*
_output_shapes
:q
SaveV2_1/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:�
SaveV2_1SaveV2ShardedFilename_1:filename:0SaveV2_1/tensor_names:output:0"SaveV2_1/shape_and_slices:output:0savev2_1_const^SaveV2"/device:CPU:0*
_output_shapes
 *
dtypes
2�
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0ShardedFilename_1:filename:0^SaveV2	^SaveV2_1"/device:CPU:0*
N*
_output_shapes
:*
T0�
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix	^SaveV2_1"/device:CPU:0*
_output_shapes
 f
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: s

Identity_1IdentityIdentity:output:0^MergeV2Checkpoints^SaveV2	^SaveV2_1*
_output_shapes
: *
T0"!

identity_1Identity_1:output:0*c
_input_shapesR
P: ::::::::: : : : : : : 2
SaveV2_1SaveV2_12
SaveV2SaveV22(
MergeV2CheckpointsMergeV2Checkpoints:	 :
 : : : : : :+ '
%
_user_specified_namefile_prefix: : : : : : : : 
�
�
F__inference_dense_600_layer_call_and_return_conditional_losses_1123777

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*'
_output_shapes
:���������*
T0N
EluEluBiasAdd:output:0*
T0*'
_output_shapes
:���������N
	Greater/yConst*
valueB
 *    *
dtype0*
_output_shapes
: j
GreaterGreaterBiasAdd:output:0Greater/y:output:0*
T0*'
_output_shapes
:���������J
mul/xConst*
valueB
 *}-�?*
dtype0*
_output_shapes
: _
mulMulmul/x:output:0Elu:activations:0*'
_output_shapes
:���������*
T0k
SelectSelectGreater:z:0Elu:activations:0mul:z:0*
T0*'
_output_shapes
:���������L
mul_1/xConst*
valueB
 *_}�?*
dtype0*
_output_shapes
: a
mul_1Mulmul_1/x:output:0Select:output:0*
T0*'
_output_shapes
:����������
IdentityIdentity	mul_1:z:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*.
_input_shapes
:���������::2.
MatMul/ReadVariableOpMatMul/ReadVariableOp20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp:& "
 
_user_specified_nameinputs: : 
�
�
+__inference_dense_601_layer_call_fn_1123809

inputs"
statefulpartitionedcall_args_1"
statefulpartitionedcall_args_2
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsstatefulpartitionedcall_args_1statefulpartitionedcall_args_2**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:���������*.
_gradient_op_typePartitionedCall-1123454*O
fJRH
F__inference_dense_601_layer_call_and_return_conditional_losses_1123448*
Tout
2�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*'
_output_shapes
:���������*
T0"
identityIdentity:output:0*.
_input_shapes
:���������::22
StatefulPartitionedCallStatefulPartitionedCall:& "
 
_user_specified_nameinputs: : 
�
�
F__inference_dense_603_layer_call_and_return_conditional_losses_1123844

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
IdentityIdentityBiasAdd:output:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:& "
 
_user_specified_nameinputs: : 
�
�
+__inference_dense_602_layer_call_fn_1123834

inputs"
statefulpartitionedcall_args_1"
statefulpartitionedcall_args_2
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsstatefulpartitionedcall_args_1statefulpartitionedcall_args_2**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:���������*.
_gradient_op_typePartitionedCall-1123489*O
fJRH
F__inference_dense_602_layer_call_and_return_conditional_losses_1123483*
Tout
2�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*.
_input_shapes
:���������::22
StatefulPartitionedCallStatefulPartitionedCall:& "
 
_user_specified_nameinputs: : 
�8
�
#__inference__traced_restore_1123973
file_prefix%
!assignvariableop_dense_600_kernel%
!assignvariableop_1_dense_600_bias'
#assignvariableop_2_dense_601_kernel%
!assignvariableop_3_dense_601_bias'
#assignvariableop_4_dense_602_kernel%
!assignvariableop_5_dense_602_bias'
#assignvariableop_6_dense_603_kernel%
!assignvariableop_7_dense_603_bias
assignvariableop_8_sgd_iter 
assignvariableop_9_sgd_decay)
%assignvariableop_10_sgd_learning_rate$
 assignvariableop_11_sgd_momentum
assignvariableop_12_total
assignvariableop_13_count
identity_15��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_10�AssignVariableOp_11�AssignVariableOp_12�AssignVariableOp_13�AssignVariableOp_2�AssignVariableOp_3�AssignVariableOp_4�AssignVariableOp_5�AssignVariableOp_6�AssignVariableOp_7�AssignVariableOp_8�AssignVariableOp_9�	RestoreV2�RestoreV2_1�
RestoreV2/tensor_namesConst"/device:CPU:0*�
value�B�B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB-optimizer/momentum/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:�
RestoreV2/shape_and_slicesConst"/device:CPU:0*/
value&B$B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:�
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*L
_output_shapes:
8::::::::::::::*
dtypes
2	L
IdentityIdentityRestoreV2:tensors:0*
_output_shapes
:*
T0}
AssignVariableOpAssignVariableOp!assignvariableop_dense_600_kernelIdentity:output:0*
dtype0*
_output_shapes
 N

Identity_1IdentityRestoreV2:tensors:1*
T0*
_output_shapes
:�
AssignVariableOp_1AssignVariableOp!assignvariableop_1_dense_600_biasIdentity_1:output:0*
_output_shapes
 *
dtype0N

Identity_2IdentityRestoreV2:tensors:2*
T0*
_output_shapes
:�
AssignVariableOp_2AssignVariableOp#assignvariableop_2_dense_601_kernelIdentity_2:output:0*
dtype0*
_output_shapes
 N

Identity_3IdentityRestoreV2:tensors:3*
T0*
_output_shapes
:�
AssignVariableOp_3AssignVariableOp!assignvariableop_3_dense_601_biasIdentity_3:output:0*
dtype0*
_output_shapes
 N

Identity_4IdentityRestoreV2:tensors:4*
T0*
_output_shapes
:�
AssignVariableOp_4AssignVariableOp#assignvariableop_4_dense_602_kernelIdentity_4:output:0*
dtype0*
_output_shapes
 N

Identity_5IdentityRestoreV2:tensors:5*
_output_shapes
:*
T0�
AssignVariableOp_5AssignVariableOp!assignvariableop_5_dense_602_biasIdentity_5:output:0*
dtype0*
_output_shapes
 N

Identity_6IdentityRestoreV2:tensors:6*
T0*
_output_shapes
:�
AssignVariableOp_6AssignVariableOp#assignvariableop_6_dense_603_kernelIdentity_6:output:0*
dtype0*
_output_shapes
 N

Identity_7IdentityRestoreV2:tensors:7*
_output_shapes
:*
T0�
AssignVariableOp_7AssignVariableOp!assignvariableop_7_dense_603_biasIdentity_7:output:0*
dtype0*
_output_shapes
 N

Identity_8IdentityRestoreV2:tensors:8*
T0	*
_output_shapes
:{
AssignVariableOp_8AssignVariableOpassignvariableop_8_sgd_iterIdentity_8:output:0*
dtype0	*
_output_shapes
 N

Identity_9IdentityRestoreV2:tensors:9*
_output_shapes
:*
T0|
AssignVariableOp_9AssignVariableOpassignvariableop_9_sgd_decayIdentity_9:output:0*
_output_shapes
 *
dtype0P
Identity_10IdentityRestoreV2:tensors:10*
T0*
_output_shapes
:�
AssignVariableOp_10AssignVariableOp%assignvariableop_10_sgd_learning_rateIdentity_10:output:0*
dtype0*
_output_shapes
 P
Identity_11IdentityRestoreV2:tensors:11*
_output_shapes
:*
T0�
AssignVariableOp_11AssignVariableOp assignvariableop_11_sgd_momentumIdentity_11:output:0*
dtype0*
_output_shapes
 P
Identity_12IdentityRestoreV2:tensors:12*
_output_shapes
:*
T0{
AssignVariableOp_12AssignVariableOpassignvariableop_12_totalIdentity_12:output:0*
dtype0*
_output_shapes
 P
Identity_13IdentityRestoreV2:tensors:13*
T0*
_output_shapes
:{
AssignVariableOp_13AssignVariableOpassignvariableop_13_countIdentity_13:output:0*
dtype0*
_output_shapes
 �
RestoreV2_1/tensor_namesConst"/device:CPU:0*1
value(B&B_CHECKPOINTABLE_OBJECT_GRAPH*
dtype0*
_output_shapes
:t
RestoreV2_1/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
valueB
B *
dtype0�
RestoreV2_1	RestoreV2file_prefix!RestoreV2_1/tensor_names:output:0%RestoreV2_1/shape_and_slices:output:0
^RestoreV2"/device:CPU:0*
_output_shapes
:*
dtypes
21
NoOpNoOp"/device:CPU:0*
_output_shapes
 �
Identity_14Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_2^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: �
Identity_15IdentityIdentity_14:output:0^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_2^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9
^RestoreV2^RestoreV2_1*
T0*
_output_shapes
: "#
identity_15Identity_15:output:0*M
_input_shapes<
:: ::::::::::::::2(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_92
	RestoreV2	RestoreV22*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122
RestoreV2_1RestoreV2_12*
AssignVariableOp_13AssignVariableOp_132(
AssignVariableOp_1AssignVariableOp_12(
AssignVariableOp_2AssignVariableOp_22(
AssignVariableOp_3AssignVariableOp_32(
AssignVariableOp_4AssignVariableOp_42$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_7: : : : : :	 :
 : : : : :+ '
%
_user_specified_namefile_prefix: : : 
�
�
F__inference_dense_602_layer_call_and_return_conditional_losses_1123483

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������N
EluEluBiasAdd:output:0*
T0*'
_output_shapes
:���������N
	Greater/yConst*
dtype0*
_output_shapes
: *
valueB
 *    j
GreaterGreaterBiasAdd:output:0Greater/y:output:0*
T0*'
_output_shapes
:���������J
mul/xConst*
valueB
 *}-�?*
dtype0*
_output_shapes
: _
mulMulmul/x:output:0Elu:activations:0*
T0*'
_output_shapes
:���������k
SelectSelectGreater:z:0Elu:activations:0mul:z:0*
T0*'
_output_shapes
:���������L
mul_1/xConst*
valueB
 *_}�?*
dtype0*
_output_shapes
: a
mul_1Mulmul_1/x:output:0Select:output:0*'
_output_shapes
:���������*
T0�
IdentityIdentity	mul_1:z:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*.
_input_shapes
:���������::2.
MatMul/ReadVariableOpMatMul/ReadVariableOp20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp:& "
 
_user_specified_nameinputs: : 
�

�
/__inference_sequential_87_layer_call_fn_1123759

inputs"
statefulpartitionedcall_args_1"
statefulpartitionedcall_args_2"
statefulpartitionedcall_args_3"
statefulpartitionedcall_args_4"
statefulpartitionedcall_args_5"
statefulpartitionedcall_args_6"
statefulpartitionedcall_args_7"
statefulpartitionedcall_args_8
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsstatefulpartitionedcall_args_1statefulpartitionedcall_args_2statefulpartitionedcall_args_3statefulpartitionedcall_args_4statefulpartitionedcall_args_5statefulpartitionedcall_args_6statefulpartitionedcall_args_7statefulpartitionedcall_args_8*.
_gradient_op_typePartitionedCall-1123598*S
fNRL
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123597*
Tout
2**
config_proto

CPU

GPU 2J 8*'
_output_shapes
:���������*
Tin
2	�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::22
StatefulPartitionedCallStatefulPartitionedCall:& "
 
_user_specified_nameinputs: : : : : : : : 
�
�
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123597

inputs,
(dense_600_statefulpartitionedcall_args_1,
(dense_600_statefulpartitionedcall_args_2,
(dense_601_statefulpartitionedcall_args_1,
(dense_601_statefulpartitionedcall_args_2,
(dense_602_statefulpartitionedcall_args_1,
(dense_602_statefulpartitionedcall_args_2,
(dense_603_statefulpartitionedcall_args_1,
(dense_603_statefulpartitionedcall_args_2
identity��!dense_600/StatefulPartitionedCall�!dense_601/StatefulPartitionedCall�!dense_602/StatefulPartitionedCall�!dense_603/StatefulPartitionedCall�
!dense_600/StatefulPartitionedCallStatefulPartitionedCallinputs(dense_600_statefulpartitionedcall_args_1(dense_600_statefulpartitionedcall_args_2*.
_gradient_op_typePartitionedCall-1123419*O
fJRH
F__inference_dense_600_layer_call_and_return_conditional_losses_1123413*
Tout
2**
config_proto

CPU

GPU 2J 8*'
_output_shapes
:���������*
Tin
2�
!dense_601/StatefulPartitionedCallStatefulPartitionedCall*dense_600/StatefulPartitionedCall:output:0(dense_601_statefulpartitionedcall_args_1(dense_601_statefulpartitionedcall_args_2*.
_gradient_op_typePartitionedCall-1123454*O
fJRH
F__inference_dense_601_layer_call_and_return_conditional_losses_1123448*
Tout
2**
config_proto

CPU

GPU 2J 8*'
_output_shapes
:���������*
Tin
2�
!dense_602/StatefulPartitionedCallStatefulPartitionedCall*dense_601/StatefulPartitionedCall:output:0(dense_602_statefulpartitionedcall_args_1(dense_602_statefulpartitionedcall_args_2*O
fJRH
F__inference_dense_602_layer_call_and_return_conditional_losses_1123483*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:���������*.
_gradient_op_typePartitionedCall-1123489�
!dense_603/StatefulPartitionedCallStatefulPartitionedCall*dense_602/StatefulPartitionedCall:output:0(dense_603_statefulpartitionedcall_args_1(dense_603_statefulpartitionedcall_args_2**
config_proto

CPU

GPU 2J 8*'
_output_shapes
:���������*
Tin
2*.
_gradient_op_typePartitionedCall-1123516*O
fJRH
F__inference_dense_603_layer_call_and_return_conditional_losses_1123510*
Tout
2�
IdentityIdentity*dense_603/StatefulPartitionedCall:output:0"^dense_600/StatefulPartitionedCall"^dense_601/StatefulPartitionedCall"^dense_602/StatefulPartitionedCall"^dense_603/StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::2F
!dense_600/StatefulPartitionedCall!dense_600/StatefulPartitionedCall2F
!dense_601/StatefulPartitionedCall!dense_601/StatefulPartitionedCall2F
!dense_602/StatefulPartitionedCall!dense_602/StatefulPartitionedCall2F
!dense_603/StatefulPartitionedCall!dense_603/StatefulPartitionedCall:& "
 
_user_specified_nameinputs: : : : : : : : 
�
�
F__inference_dense_601_layer_call_and_return_conditional_losses_1123802

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*'
_output_shapes
:���������*
T0�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������N
EluEluBiasAdd:output:0*
T0*'
_output_shapes
:���������N
	Greater/yConst*
valueB
 *    *
dtype0*
_output_shapes
: j
GreaterGreaterBiasAdd:output:0Greater/y:output:0*
T0*'
_output_shapes
:���������J
mul/xConst*
_output_shapes
: *
valueB
 *}-�?*
dtype0_
mulMulmul/x:output:0Elu:activations:0*'
_output_shapes
:���������*
T0k
SelectSelectGreater:z:0Elu:activations:0mul:z:0*
T0*'
_output_shapes
:���������L
mul_1/xConst*
valueB
 *_}�?*
dtype0*
_output_shapes
: a
mul_1Mulmul_1/x:output:0Select:output:0*
T0*'
_output_shapes
:����������
IdentityIdentity	mul_1:z:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*'
_output_shapes
:���������*
T0"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:& "
 
_user_specified_nameinputs: : 
�
�
F__inference_dense_602_layer_call_and_return_conditional_losses_1123827

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������N
EluEluBiasAdd:output:0*
T0*'
_output_shapes
:���������N
	Greater/yConst*
_output_shapes
: *
valueB
 *    *
dtype0j
GreaterGreaterBiasAdd:output:0Greater/y:output:0*
T0*'
_output_shapes
:���������J
mul/xConst*
_output_shapes
: *
valueB
 *}-�?*
dtype0_
mulMulmul/x:output:0Elu:activations:0*
T0*'
_output_shapes
:���������k
SelectSelectGreater:z:0Elu:activations:0mul:z:0*
T0*'
_output_shapes
:���������L
mul_1/xConst*
valueB
 *_}�?*
dtype0*
_output_shapes
: a
mul_1Mulmul_1/x:output:0Select:output:0*
T0*'
_output_shapes
:����������
IdentityIdentity	mul_1:z:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*.
_input_shapes
:���������::2.
MatMul/ReadVariableOpMatMul/ReadVariableOp20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp:& "
 
_user_specified_nameinputs: : 
�
�
F__inference_dense_603_layer_call_and_return_conditional_losses_1123510

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
IdentityIdentityBiasAdd:output:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*.
_input_shapes
:���������::2.
MatMul/ReadVariableOpMatMul/ReadVariableOp20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp: : :& "
 
_user_specified_nameinputs
�
�
F__inference_dense_601_layer_call_and_return_conditional_losses_1123448

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������N
EluEluBiasAdd:output:0*'
_output_shapes
:���������*
T0N
	Greater/yConst*
valueB
 *    *
dtype0*
_output_shapes
: j
GreaterGreaterBiasAdd:output:0Greater/y:output:0*
T0*'
_output_shapes
:���������J
mul/xConst*
valueB
 *}-�?*
dtype0*
_output_shapes
: _
mulMulmul/x:output:0Elu:activations:0*
T0*'
_output_shapes
:���������k
SelectSelectGreater:z:0Elu:activations:0mul:z:0*
T0*'
_output_shapes
:���������L
mul_1/xConst*
valueB
 *_}�?*
dtype0*
_output_shapes
: a
mul_1Mulmul_1/x:output:0Select:output:0*
T0*'
_output_shapes
:����������
IdentityIdentity	mul_1:z:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp: :& "
 
_user_specified_nameinputs: 
�
�
+__inference_dense_600_layer_call_fn_1123784

inputs"
statefulpartitionedcall_args_1"
statefulpartitionedcall_args_2
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsstatefulpartitionedcall_args_1statefulpartitionedcall_args_2*.
_gradient_op_typePartitionedCall-1123419*O
fJRH
F__inference_dense_600_layer_call_and_return_conditional_losses_1123413*
Tout
2**
config_proto

CPU

GPU 2J 8*'
_output_shapes
:���������*
Tin
2�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*'
_output_shapes
:���������*
T0"
identityIdentity:output:0*.
_input_shapes
:���������::22
StatefulPartitionedCallStatefulPartitionedCall:& "
 
_user_specified_nameinputs: : 
�
�
+__inference_dense_603_layer_call_fn_1123851

inputs"
statefulpartitionedcall_args_1"
statefulpartitionedcall_args_2
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsstatefulpartitionedcall_args_1statefulpartitionedcall_args_2*
Tout
2**
config_proto

CPU

GPU 2J 8*'
_output_shapes
:���������*
Tin
2*.
_gradient_op_typePartitionedCall-1123516*O
fJRH
F__inference_dense_603_layer_call_and_return_conditional_losses_1123510�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*'
_output_shapes
:���������*
T0"
identityIdentity:output:0*.
_input_shapes
:���������::22
StatefulPartitionedCallStatefulPartitionedCall:& "
 
_user_specified_nameinputs: : 
�8
�
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123681

inputs,
(dense_600_matmul_readvariableop_resource-
)dense_600_biasadd_readvariableop_resource,
(dense_601_matmul_readvariableop_resource-
)dense_601_biasadd_readvariableop_resource,
(dense_602_matmul_readvariableop_resource-
)dense_602_biasadd_readvariableop_resource,
(dense_603_matmul_readvariableop_resource-
)dense_603_biasadd_readvariableop_resource
identity�� dense_600/BiasAdd/ReadVariableOp�dense_600/MatMul/ReadVariableOp� dense_601/BiasAdd/ReadVariableOp�dense_601/MatMul/ReadVariableOp� dense_602/BiasAdd/ReadVariableOp�dense_602/MatMul/ReadVariableOp� dense_603/BiasAdd/ReadVariableOp�dense_603/MatMul/ReadVariableOp�
dense_600/MatMul/ReadVariableOpReadVariableOp(dense_600_matmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:}
dense_600/MatMulMatMulinputs'dense_600/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
 dense_600/BiasAdd/ReadVariableOpReadVariableOp)dense_600_biasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:�
dense_600/BiasAddBiasAdddense_600/MatMul:product:0(dense_600/BiasAdd/ReadVariableOp:value:0*'
_output_shapes
:���������*
T0b
dense_600/EluEludense_600/BiasAdd:output:0*
T0*'
_output_shapes
:���������X
dense_600/Greater/yConst*
valueB
 *    *
dtype0*
_output_shapes
: �
dense_600/GreaterGreaterdense_600/BiasAdd:output:0dense_600/Greater/y:output:0*'
_output_shapes
:���������*
T0T
dense_600/mul/xConst*
_output_shapes
: *
valueB
 *}-�?*
dtype0}
dense_600/mulMuldense_600/mul/x:output:0dense_600/Elu:activations:0*
T0*'
_output_shapes
:����������
dense_600/SelectSelectdense_600/Greater:z:0dense_600/Elu:activations:0dense_600/mul:z:0*'
_output_shapes
:���������*
T0V
dense_600/mul_1/xConst*
_output_shapes
: *
valueB
 *_}�?*
dtype0
dense_600/mul_1Muldense_600/mul_1/x:output:0dense_600/Select:output:0*
T0*'
_output_shapes
:����������
dense_601/MatMul/ReadVariableOpReadVariableOp(dense_601_matmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:�
dense_601/MatMulMatMuldense_600/mul_1:z:0'dense_601/MatMul/ReadVariableOp:value:0*'
_output_shapes
:���������*
T0�
 dense_601/BiasAdd/ReadVariableOpReadVariableOp)dense_601_biasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:�
dense_601/BiasAddBiasAdddense_601/MatMul:product:0(dense_601/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������b
dense_601/EluEludense_601/BiasAdd:output:0*'
_output_shapes
:���������*
T0X
dense_601/Greater/yConst*
valueB
 *    *
dtype0*
_output_shapes
: �
dense_601/GreaterGreaterdense_601/BiasAdd:output:0dense_601/Greater/y:output:0*
T0*'
_output_shapes
:���������T
dense_601/mul/xConst*
valueB
 *}-�?*
dtype0*
_output_shapes
: }
dense_601/mulMuldense_601/mul/x:output:0dense_601/Elu:activations:0*
T0*'
_output_shapes
:����������
dense_601/SelectSelectdense_601/Greater:z:0dense_601/Elu:activations:0dense_601/mul:z:0*
T0*'
_output_shapes
:���������V
dense_601/mul_1/xConst*
_output_shapes
: *
valueB
 *_}�?*
dtype0
dense_601/mul_1Muldense_601/mul_1/x:output:0dense_601/Select:output:0*
T0*'
_output_shapes
:����������
dense_602/MatMul/ReadVariableOpReadVariableOp(dense_602_matmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:�
dense_602/MatMulMatMuldense_601/mul_1:z:0'dense_602/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
 dense_602/BiasAdd/ReadVariableOpReadVariableOp)dense_602_biasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:�
dense_602/BiasAddBiasAdddense_602/MatMul:product:0(dense_602/BiasAdd/ReadVariableOp:value:0*'
_output_shapes
:���������*
T0b
dense_602/EluEludense_602/BiasAdd:output:0*
T0*'
_output_shapes
:���������X
dense_602/Greater/yConst*
_output_shapes
: *
valueB
 *    *
dtype0�
dense_602/GreaterGreaterdense_602/BiasAdd:output:0dense_602/Greater/y:output:0*'
_output_shapes
:���������*
T0T
dense_602/mul/xConst*
valueB
 *}-�?*
dtype0*
_output_shapes
: }
dense_602/mulMuldense_602/mul/x:output:0dense_602/Elu:activations:0*'
_output_shapes
:���������*
T0�
dense_602/SelectSelectdense_602/Greater:z:0dense_602/Elu:activations:0dense_602/mul:z:0*'
_output_shapes
:���������*
T0V
dense_602/mul_1/xConst*
valueB
 *_}�?*
dtype0*
_output_shapes
: 
dense_602/mul_1Muldense_602/mul_1/x:output:0dense_602/Select:output:0*'
_output_shapes
:���������*
T0�
dense_603/MatMul/ReadVariableOpReadVariableOp(dense_603_matmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:�
dense_603/MatMulMatMuldense_602/mul_1:z:0'dense_603/MatMul/ReadVariableOp:value:0*'
_output_shapes
:���������*
T0�
 dense_603/BiasAdd/ReadVariableOpReadVariableOp)dense_603_biasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:�
dense_603/BiasAddBiasAdddense_603/MatMul:product:0(dense_603/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
IdentityIdentitydense_603/BiasAdd:output:0!^dense_600/BiasAdd/ReadVariableOp ^dense_600/MatMul/ReadVariableOp!^dense_601/BiasAdd/ReadVariableOp ^dense_601/MatMul/ReadVariableOp!^dense_602/BiasAdd/ReadVariableOp ^dense_602/MatMul/ReadVariableOp!^dense_603/BiasAdd/ReadVariableOp ^dense_603/MatMul/ReadVariableOp*'
_output_shapes
:���������*
T0"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::2B
dense_603/MatMul/ReadVariableOpdense_603/MatMul/ReadVariableOp2D
 dense_603/BiasAdd/ReadVariableOp dense_603/BiasAdd/ReadVariableOp2D
 dense_602/BiasAdd/ReadVariableOp dense_602/BiasAdd/ReadVariableOp2B
dense_600/MatMul/ReadVariableOpdense_600/MatMul/ReadVariableOp2D
 dense_601/BiasAdd/ReadVariableOp dense_601/BiasAdd/ReadVariableOp2B
dense_602/MatMul/ReadVariableOpdense_602/MatMul/ReadVariableOp2D
 dense_600/BiasAdd/ReadVariableOp dense_600/BiasAdd/ReadVariableOp2B
dense_601/MatMul/ReadVariableOpdense_601/MatMul/ReadVariableOp: : :& "
 
_user_specified_nameinputs: : : : : : 
�G
�
"__inference__wrapped_model_1123389
dense_600_input:
6sequential_87_dense_600_matmul_readvariableop_resource;
7sequential_87_dense_600_biasadd_readvariableop_resource:
6sequential_87_dense_601_matmul_readvariableop_resource;
7sequential_87_dense_601_biasadd_readvariableop_resource:
6sequential_87_dense_602_matmul_readvariableop_resource;
7sequential_87_dense_602_biasadd_readvariableop_resource:
6sequential_87_dense_603_matmul_readvariableop_resource;
7sequential_87_dense_603_biasadd_readvariableop_resource
identity��.sequential_87/dense_600/BiasAdd/ReadVariableOp�-sequential_87/dense_600/MatMul/ReadVariableOp�.sequential_87/dense_601/BiasAdd/ReadVariableOp�-sequential_87/dense_601/MatMul/ReadVariableOp�.sequential_87/dense_602/BiasAdd/ReadVariableOp�-sequential_87/dense_602/MatMul/ReadVariableOp�.sequential_87/dense_603/BiasAdd/ReadVariableOp�-sequential_87/dense_603/MatMul/ReadVariableOp�
-sequential_87/dense_600/MatMul/ReadVariableOpReadVariableOp6sequential_87_dense_600_matmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:�
sequential_87/dense_600/MatMulMatMuldense_600_input5sequential_87/dense_600/MatMul/ReadVariableOp:value:0*'
_output_shapes
:���������*
T0�
.sequential_87/dense_600/BiasAdd/ReadVariableOpReadVariableOp7sequential_87_dense_600_biasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:�
sequential_87/dense_600/BiasAddBiasAdd(sequential_87/dense_600/MatMul:product:06sequential_87/dense_600/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������~
sequential_87/dense_600/EluElu(sequential_87/dense_600/BiasAdd:output:0*
T0*'
_output_shapes
:���������f
!sequential_87/dense_600/Greater/yConst*
valueB
 *    *
dtype0*
_output_shapes
: �
sequential_87/dense_600/GreaterGreater(sequential_87/dense_600/BiasAdd:output:0*sequential_87/dense_600/Greater/y:output:0*
T0*'
_output_shapes
:���������b
sequential_87/dense_600/mul/xConst*
dtype0*
_output_shapes
: *
valueB
 *}-�?�
sequential_87/dense_600/mulMul&sequential_87/dense_600/mul/x:output:0)sequential_87/dense_600/Elu:activations:0*'
_output_shapes
:���������*
T0�
sequential_87/dense_600/SelectSelect#sequential_87/dense_600/Greater:z:0)sequential_87/dense_600/Elu:activations:0sequential_87/dense_600/mul:z:0*
T0*'
_output_shapes
:���������d
sequential_87/dense_600/mul_1/xConst*
valueB
 *_}�?*
dtype0*
_output_shapes
: �
sequential_87/dense_600/mul_1Mul(sequential_87/dense_600/mul_1/x:output:0'sequential_87/dense_600/Select:output:0*
T0*'
_output_shapes
:����������
-sequential_87/dense_601/MatMul/ReadVariableOpReadVariableOp6sequential_87_dense_601_matmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:�
sequential_87/dense_601/MatMulMatMul!sequential_87/dense_600/mul_1:z:05sequential_87/dense_601/MatMul/ReadVariableOp:value:0*'
_output_shapes
:���������*
T0�
.sequential_87/dense_601/BiasAdd/ReadVariableOpReadVariableOp7sequential_87_dense_601_biasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:�
sequential_87/dense_601/BiasAddBiasAdd(sequential_87/dense_601/MatMul:product:06sequential_87/dense_601/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������~
sequential_87/dense_601/EluElu(sequential_87/dense_601/BiasAdd:output:0*
T0*'
_output_shapes
:���������f
!sequential_87/dense_601/Greater/yConst*
valueB
 *    *
dtype0*
_output_shapes
: �
sequential_87/dense_601/GreaterGreater(sequential_87/dense_601/BiasAdd:output:0*sequential_87/dense_601/Greater/y:output:0*
T0*'
_output_shapes
:���������b
sequential_87/dense_601/mul/xConst*
valueB
 *}-�?*
dtype0*
_output_shapes
: �
sequential_87/dense_601/mulMul&sequential_87/dense_601/mul/x:output:0)sequential_87/dense_601/Elu:activations:0*
T0*'
_output_shapes
:����������
sequential_87/dense_601/SelectSelect#sequential_87/dense_601/Greater:z:0)sequential_87/dense_601/Elu:activations:0sequential_87/dense_601/mul:z:0*
T0*'
_output_shapes
:���������d
sequential_87/dense_601/mul_1/xConst*
valueB
 *_}�?*
dtype0*
_output_shapes
: �
sequential_87/dense_601/mul_1Mul(sequential_87/dense_601/mul_1/x:output:0'sequential_87/dense_601/Select:output:0*'
_output_shapes
:���������*
T0�
-sequential_87/dense_602/MatMul/ReadVariableOpReadVariableOp6sequential_87_dense_602_matmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:�
sequential_87/dense_602/MatMulMatMul!sequential_87/dense_601/mul_1:z:05sequential_87/dense_602/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
.sequential_87/dense_602/BiasAdd/ReadVariableOpReadVariableOp7sequential_87_dense_602_biasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:�
sequential_87/dense_602/BiasAddBiasAdd(sequential_87/dense_602/MatMul:product:06sequential_87/dense_602/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������~
sequential_87/dense_602/EluElu(sequential_87/dense_602/BiasAdd:output:0*'
_output_shapes
:���������*
T0f
!sequential_87/dense_602/Greater/yConst*
valueB
 *    *
dtype0*
_output_shapes
: �
sequential_87/dense_602/GreaterGreater(sequential_87/dense_602/BiasAdd:output:0*sequential_87/dense_602/Greater/y:output:0*'
_output_shapes
:���������*
T0b
sequential_87/dense_602/mul/xConst*
valueB
 *}-�?*
dtype0*
_output_shapes
: �
sequential_87/dense_602/mulMul&sequential_87/dense_602/mul/x:output:0)sequential_87/dense_602/Elu:activations:0*
T0*'
_output_shapes
:����������
sequential_87/dense_602/SelectSelect#sequential_87/dense_602/Greater:z:0)sequential_87/dense_602/Elu:activations:0sequential_87/dense_602/mul:z:0*'
_output_shapes
:���������*
T0d
sequential_87/dense_602/mul_1/xConst*
valueB
 *_}�?*
dtype0*
_output_shapes
: �
sequential_87/dense_602/mul_1Mul(sequential_87/dense_602/mul_1/x:output:0'sequential_87/dense_602/Select:output:0*
T0*'
_output_shapes
:����������
-sequential_87/dense_603/MatMul/ReadVariableOpReadVariableOp6sequential_87_dense_603_matmul_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes

:�
sequential_87/dense_603/MatMulMatMul!sequential_87/dense_602/mul_1:z:05sequential_87/dense_603/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
.sequential_87/dense_603/BiasAdd/ReadVariableOpReadVariableOp7sequential_87_dense_603_biasadd_readvariableop_resource",/job:localhost/replica:0/task:0/device:CPU:0*
dtype0*
_output_shapes
:�
sequential_87/dense_603/BiasAddBiasAdd(sequential_87/dense_603/MatMul:product:06sequential_87/dense_603/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
IdentityIdentity(sequential_87/dense_603/BiasAdd:output:0/^sequential_87/dense_600/BiasAdd/ReadVariableOp.^sequential_87/dense_600/MatMul/ReadVariableOp/^sequential_87/dense_601/BiasAdd/ReadVariableOp.^sequential_87/dense_601/MatMul/ReadVariableOp/^sequential_87/dense_602/BiasAdd/ReadVariableOp.^sequential_87/dense_602/MatMul/ReadVariableOp/^sequential_87/dense_603/BiasAdd/ReadVariableOp.^sequential_87/dense_603/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::2^
-sequential_87/dense_602/MatMul/ReadVariableOp-sequential_87/dense_602/MatMul/ReadVariableOp2`
.sequential_87/dense_600/BiasAdd/ReadVariableOp.sequential_87/dense_600/BiasAdd/ReadVariableOp2^
-sequential_87/dense_601/MatMul/ReadVariableOp-sequential_87/dense_601/MatMul/ReadVariableOp2`
.sequential_87/dense_603/BiasAdd/ReadVariableOp.sequential_87/dense_603/BiasAdd/ReadVariableOp2^
-sequential_87/dense_603/MatMul/ReadVariableOp-sequential_87/dense_603/MatMul/ReadVariableOp2`
.sequential_87/dense_602/BiasAdd/ReadVariableOp.sequential_87/dense_602/BiasAdd/ReadVariableOp2^
-sequential_87/dense_600/MatMul/ReadVariableOp-sequential_87/dense_600/MatMul/ReadVariableOp2`
.sequential_87/dense_601/BiasAdd/ReadVariableOp.sequential_87/dense_601/BiasAdd/ReadVariableOp: : :/ +
)
_user_specified_namedense_600_input: : : : : : 
�

�
%__inference_signature_wrapper_1123627
dense_600_input"
statefulpartitionedcall_args_1"
statefulpartitionedcall_args_2"
statefulpartitionedcall_args_3"
statefulpartitionedcall_args_4"
statefulpartitionedcall_args_5"
statefulpartitionedcall_args_6"
statefulpartitionedcall_args_7"
statefulpartitionedcall_args_8
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCalldense_600_inputstatefulpartitionedcall_args_1statefulpartitionedcall_args_2statefulpartitionedcall_args_3statefulpartitionedcall_args_4statefulpartitionedcall_args_5statefulpartitionedcall_args_6statefulpartitionedcall_args_7statefulpartitionedcall_args_8*.
_gradient_op_typePartitionedCall-1123616*+
f&R$
"__inference__wrapped_model_1123389*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2	*'
_output_shapes
:����������
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*F
_input_shapes5
3:���������::::::::22
StatefulPartitionedCallStatefulPartitionedCall: : : : : : :/ +
)
_user_specified_namedense_600_input: : "wL
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*�
serving_default�
K
dense_600_input8
!serving_default_dense_600_input:0���������=
	dense_6030
StatefulPartitionedCall:0���������tensorflow/serving/predict*>
__saved_model_init_op%#
__saved_model_init_op

NoOp:��
�$
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer_with_weights-3
layer-4
	optimizer
trainable_variables
	variables
	regularization_losses

	keras_api

signatures
P_default_save_signature
*Q&call_and_return_all_conditional_losses
R__call__"�!
_tf_keras_sequential�!{"class_name": "Sequential", "name": "sequential_87", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "config": {"name": "sequential_87", "layers": [{"class_name": "Dense", "config": {"name": "dense_600", "trainable": true, "batch_input_shape": [null, 1], "dtype": "float32", "units": 1, "activation": "selu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_601", "trainable": true, "dtype": "float32", "units": 1, "activation": "selu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_602", "trainable": true, "dtype": "float32", "units": 1, "activation": "selu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_603", "trainable": true, "dtype": "float32", "units": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}]}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 1}}}, "keras_version": "2.2.4-tf", "backend": "tensorflow", "model_config": {"class_name": "Sequential", "config": {"name": "sequential_87", "layers": [{"class_name": "Dense", "config": {"name": "dense_600", "trainable": true, "batch_input_shape": [null, 1], "dtype": "float32", "units": 1, "activation": "selu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_601", "trainable": true, "dtype": "float32", "units": 1, "activation": "selu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_602", "trainable": true, "dtype": "float32", "units": 1, "activation": "selu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_603", "trainable": true, "dtype": "float32", "units": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}]}}, "training_config": {"loss": "mse", "metrics": ["mse"], "weighted_metrics": null, "sample_weight_mode": null, "loss_weights": null, "optimizer_config": {"class_name": "SGD", "config": {"name": "SGD", "learning_rate": 0.009999999776482582, "decay": 0.0, "momentum": 0.0, "nesterov": false}}}}
�
trainable_variables
	variables
regularization_losses
	keras_api
*S&call_and_return_all_conditional_losses
T__call__"�
_tf_keras_layer�{"class_name": "InputLayer", "name": "dense_600_input", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": [null, 1], "config": {"batch_input_shape": [null, 1], "dtype": "float32", "sparse": false, "name": "dense_600_input"}}
�

kernel
bias
trainable_variables
	variables
regularization_losses
	keras_api
*U&call_and_return_all_conditional_losses
V__call__"�
_tf_keras_layer�{"class_name": "Dense", "name": "dense_600", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": [null, 1], "config": {"name": "dense_600", "trainable": true, "batch_input_shape": [null, 1], "dtype": "float32", "units": 1, "activation": "selu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 1}}}}
�

kernel
bias
trainable_variables
	variables
regularization_losses
	keras_api
*W&call_and_return_all_conditional_losses
X__call__"�
_tf_keras_layer�{"class_name": "Dense", "name": "dense_601", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "config": {"name": "dense_601", "trainable": true, "dtype": "float32", "units": 1, "activation": "selu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 1}}}}
�

kernel
bias
trainable_variables
	variables
 regularization_losses
!	keras_api
*Y&call_and_return_all_conditional_losses
Z__call__"�
_tf_keras_layer�{"class_name": "Dense", "name": "dense_602", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "config": {"name": "dense_602", "trainable": true, "dtype": "float32", "units": 1, "activation": "selu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 1}}}}
�

"kernel
#bias
$trainable_variables
%	variables
&regularization_losses
'	keras_api
*[&call_and_return_all_conditional_losses
\__call__"�
_tf_keras_layer�{"class_name": "Dense", "name": "dense_603", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "config": {"name": "dense_603", "trainable": true, "dtype": "float32", "units": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 1}}}}
I
(iter
	)decay
*learning_rate
+momentum"
	optimizer
X
0
1
2
3
4
5
"6
#7"
trackable_list_wrapper
X
0
1
2
3
4
5
"6
#7"
trackable_list_wrapper
 "
trackable_list_wrapper
�

,layers
-layer_regularization_losses
.non_trainable_variables
trainable_variables
/metrics
	variables
	regularization_losses
R__call__
P_default_save_signature
*Q&call_and_return_all_conditional_losses
&Q"call_and_return_conditional_losses"
_generic_user_object
,
]serving_default"
signature_map
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�

0layers
1layer_regularization_losses
2non_trainable_variables
trainable_variables
3metrics
	variables
regularization_losses
T__call__
*S&call_and_return_all_conditional_losses
&S"call_and_return_conditional_losses"
_generic_user_object
": 2dense_600/kernel
:2dense_600/bias
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�

4layers
5layer_regularization_losses
6non_trainable_variables
trainable_variables
7metrics
	variables
regularization_losses
V__call__
*U&call_and_return_all_conditional_losses
&U"call_and_return_conditional_losses"
_generic_user_object
": 2dense_601/kernel
:2dense_601/bias
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�

8layers
9layer_regularization_losses
:non_trainable_variables
trainable_variables
;metrics
	variables
regularization_losses
X__call__
*W&call_and_return_all_conditional_losses
&W"call_and_return_conditional_losses"
_generic_user_object
": 2dense_602/kernel
:2dense_602/bias
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�

<layers
=layer_regularization_losses
>non_trainable_variables
trainable_variables
?metrics
	variables
 regularization_losses
Z__call__
*Y&call_and_return_all_conditional_losses
&Y"call_and_return_conditional_losses"
_generic_user_object
": 2dense_603/kernel
:2dense_603/bias
.
"0
#1"
trackable_list_wrapper
.
"0
#1"
trackable_list_wrapper
 "
trackable_list_wrapper
�

@layers
Alayer_regularization_losses
Bnon_trainable_variables
$trainable_variables
Cmetrics
%	variables
&regularization_losses
\__call__
*[&call_and_return_all_conditional_losses
&["call_and_return_conditional_losses"
_generic_user_object
:	 (2SGD/iter
: (2	SGD/decay
: (2SGD/learning_rate
: (2SGD/momentum
<
0
1
2
3"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
'
D0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
	Etotal
	Fcount
G
_fn_kwargs
Htrainable_variables
I	variables
Jregularization_losses
K	keras_api
*^&call_and_return_all_conditional_losses
___call__"�
_tf_keras_layer�{"class_name": "MeanMetricWrapper", "name": "mse", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "config": {"name": "mse", "dtype": "float32"}}
:  (2total
:  (2count
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
.
E0
F1"
trackable_list_wrapper
 "
trackable_list_wrapper
�

Llayers
Mlayer_regularization_losses
Nnon_trainable_variables
Htrainable_variables
Ometrics
I	variables
Jregularization_losses
___call__
*^&call_and_return_all_conditional_losses
&^"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
E0
F1"
trackable_list_wrapper
 "
trackable_list_wrapper
�2�
"__inference__wrapped_model_1123389�
���
FullArgSpec
args� 
varargsjargs
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *.�+
)�&
dense_600_input���������
�2�
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123681
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123546
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123733
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123528�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
/__inference_sequential_87_layer_call_fn_1123577
/__inference_sequential_87_layer_call_fn_1123746
/__inference_sequential_87_layer_call_fn_1123609
/__inference_sequential_87_layer_call_fn_1123759�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2��
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkwjkwargs
defaults� 

kwonlyargs�

jtraining%
kwonlydefaults�

trainingp 
annotations� *
 
�2��
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkwjkwargs
defaults� 

kwonlyargs�

jtraining%
kwonlydefaults�

trainingp 
annotations� *
 
�2�
F__inference_dense_600_layer_call_and_return_conditional_losses_1123777�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_dense_600_layer_call_fn_1123784�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_dense_601_layer_call_and_return_conditional_losses_1123802�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_dense_601_layer_call_fn_1123809�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_dense_602_layer_call_and_return_conditional_losses_1123827�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_dense_602_layer_call_fn_1123834�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_dense_603_layer_call_and_return_conditional_losses_1123844�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_dense_603_layer_call_fn_1123851�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
<B:
%__inference_signature_wrapper_1123627dense_600_input
�2��
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkwjkwargs
defaults� 

kwonlyargs�

jtraining%
kwonlydefaults�

trainingp 
annotations� *
 
�2��
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkwjkwargs
defaults� 

kwonlyargs�

jtraining%
kwonlydefaults�

trainingp 
annotations� *
 �
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123681j"#7�4
-�*
 �
inputs���������
p

 
� "%�"
�
0���������
� �
F__inference_dense_600_layer_call_and_return_conditional_losses_1123777\/�,
%�"
 �
inputs���������
� "%�"
�
0���������
� �
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123733j"#7�4
-�*
 �
inputs���������
p 

 
� "%�"
�
0���������
� �
/__inference_sequential_87_layer_call_fn_1123609f"#@�=
6�3
)�&
dense_600_input���������
p 

 
� "�����������
F__inference_dense_602_layer_call_and_return_conditional_losses_1123827\/�,
%�"
 �
inputs���������
� "%�"
�
0���������
� ~
+__inference_dense_602_layer_call_fn_1123834O/�,
%�"
 �
inputs���������
� "�����������
%__inference_signature_wrapper_1123627�"#K�H
� 
A�>
<
dense_600_input)�&
dense_600_input���������"5�2
0
	dense_603#� 
	dense_603����������
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123528s"#@�=
6�3
)�&
dense_600_input���������
p

 
� "%�"
�
0���������
� �
/__inference_sequential_87_layer_call_fn_1123577f"#@�=
6�3
)�&
dense_600_input���������
p

 
� "�����������
/__inference_sequential_87_layer_call_fn_1123746]"#7�4
-�*
 �
inputs���������
p

 
� "�����������
F__inference_dense_603_layer_call_and_return_conditional_losses_1123844\"#/�,
%�"
 �
inputs���������
� "%�"
�
0���������
� �
F__inference_dense_601_layer_call_and_return_conditional_losses_1123802\/�,
%�"
 �
inputs���������
� "%�"
�
0���������
� �
J__inference_sequential_87_layer_call_and_return_conditional_losses_1123546s"#@�=
6�3
)�&
dense_600_input���������
p 

 
� "%�"
�
0���������
� ~
+__inference_dense_601_layer_call_fn_1123809O/�,
%�"
 �
inputs���������
� "����������~
+__inference_dense_600_layer_call_fn_1123784O/�,
%�"
 �
inputs���������
� "�����������
"__inference__wrapped_model_1123389{"#8�5
.�+
)�&
dense_600_input���������
� "5�2
0
	dense_603#� 
	dense_603���������~
+__inference_dense_603_layer_call_fn_1123851O"#/�,
%�"
 �
inputs���������
� "�����������
/__inference_sequential_87_layer_call_fn_1123759]"#7�4
-�*
 �
inputs���������
p 

 
� "����������