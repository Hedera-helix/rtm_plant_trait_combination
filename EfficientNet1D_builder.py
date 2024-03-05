import math


import tensorflow as tf
print(tf.__version__)

from tensorflow.keras.models import Sequential
from tensorflow.keras.models import Model


#######
from tensorflow.keras import layers
from tensorflow.keras.layers import InputLayer,Input
from tensorflow.keras.layers import Conv1D,MaxPooling1D,GlobalAveragePooling1D,AveragePooling1D,ZeroPadding1D
from tensorflow.keras.layers import Dense,Reshape,Dropout,Flatten,BatchNormalization


import string
import collections
from six.moves import xrange



BlockArgs = collections.namedtuple('BlockArgs', [
    'kernel_size', 'num_repeat', 'input_filters', 'output_filters',
    'expand_ratio', 'id_skip', 'strides', 'se_ratio'
])
# defaults will be a public argument for namedtuple in Python 3.7
# https://docs.python.org/3/library/collections.html#collections.namedtuple
BlockArgs.__new__.__defaults__ = (None,) * len(BlockArgs._fields)

DEFAULT_BLOCKS_ARGS = [
    BlockArgs(kernel_size=3, num_repeat=1, input_filters=32, output_filters=16,
              expand_ratio=1, id_skip=True, strides=[1], se_ratio=0.25),
    BlockArgs(kernel_size=3, num_repeat=2, input_filters=16, output_filters=24,
              expand_ratio=6, id_skip=True, strides=[2], se_ratio=0.25),
    BlockArgs(kernel_size=5, num_repeat=2, input_filters=24, output_filters=40,
              expand_ratio=6, id_skip=True, strides=[2], se_ratio=0.25),
    BlockArgs(kernel_size=3, num_repeat=3, input_filters=40, output_filters=80,
              expand_ratio=6, id_skip=True, strides=[2], se_ratio=0.25),
    BlockArgs(kernel_size=5, num_repeat=3, input_filters=80, output_filters=112,
              expand_ratio=6, id_skip=True, strides=[1], se_ratio=0.25),
    BlockArgs(kernel_size=5, num_repeat=4, input_filters=112, output_filters=192,
              expand_ratio=6, id_skip=True, strides=[2], se_ratio=0.25),
    BlockArgs(kernel_size=3, num_repeat=1, input_filters=192, output_filters=320,
              expand_ratio=6, id_skip=True, strides=[1], se_ratio=0.25)
]
CONV_KERNEL_INITIALIZER = {
    'class_name': 'VarianceScaling',
    'config': {
        'scale': 2.0,
        'mode': 'fan_out',
        # EfficientNet actually uses an untruncated normal distribution for
        # initializing conv layers, but keras.initializers.VarianceScaling use
        # a truncated distribution.
        # We decided against a custom initializer for better serializability.
        'distribution': 'normal'
    }
}


def mb_conv_block_1d(inputs, block_args, activation, drop_rate=None, prefix='', ):
    """Mobile Inverted Residual Bottleneck."""

    has_se = (block_args.se_ratio is not None) and (0 < block_args.se_ratio <= 1)
    bn_axis = 2

    # Expansion phase
    filters = block_args.input_filters * block_args.expand_ratio
    if block_args.expand_ratio != 1:
        x = layers.Conv1D(filters, 1,
                          padding='same',
                          use_bias=False,
                          kernel_initializer=CONV_KERNEL_INITIALIZER,
                          name=prefix + 'expand_conv')(inputs)
        x = layers.BatchNormalization(axis=bn_axis, name=prefix + 'expand_bn')(x)
        x = layers.Activation(activation, name=prefix + 'expand_activation')(x)
    else:
        x = inputs

    # Depthwise Convolution
    x = layers.DepthwiseConv1D(block_args.kernel_size,
                               strides=block_args.strides,
                               padding='same',
                               use_bias=False,
                               depthwise_initializer=CONV_KERNEL_INITIALIZER,
                               name=prefix + 'dwconv')(x)
    x = layers.BatchNormalization(axis=bn_axis, name=prefix + 'bn')(x)
    x = layers.Activation(activation, name=prefix + 'activation')(x)

    # Squeeze and Excitation phase
    if has_se:
        num_reduced_filters = max(1, int(
            block_args.input_filters * block_args.se_ratio
        ))
        se_tensor = layers.GlobalAveragePooling1D(name=prefix + 'se_squeeze')(x)

        target_shape = (1, filters)
        se_tensor = layers.Reshape(target_shape, name=prefix + 'se_reshape')(se_tensor)
        se_tensor = layers.Conv1D(num_reduced_filters, 1,
                                  activation=activation,
                                  padding='same',
                                  use_bias=True,
                                  kernel_initializer=CONV_KERNEL_INITIALIZER,
                                  name=prefix + 'se_reduce')(se_tensor)
        se_tensor = layers.Conv1D(filters, 1,
                                  activation='sigmoid',
                                  padding='same',
                                  use_bias=True,
                                  kernel_initializer=CONV_KERNEL_INITIALIZER,
                                  name=prefix + 'se_expand')(se_tensor)

        x = layers.multiply([x, se_tensor], name=prefix + 'se_excite')

    # Output phase
    x = layers.Conv1D(block_args.output_filters, 1,
                      padding='same',
                      use_bias=False,
                      kernel_initializer=CONV_KERNEL_INITIALIZER,
                      name=prefix + 'project_conv')(x)
    x = layers.BatchNormalization(axis=bn_axis, name=prefix + 'project_bn')(x)
    
    if block_args.id_skip and block_args.strides == 1 and block_args.input_filters == block_args.output_filters: #and all(s == 1 for s in block_args.strides)
        if drop_rate and (drop_rate > 0):
            x = Dropout(drop_rate,name=prefix + 'drop')(x)
#             x = Dropout(drop_rate,
#                         noise_shape=(None, 1, 1, 1),
#                         name=prefix + 'drop')(x)
        x = layers.add([x, inputs], name=prefix + 'add')

    return x

def round_filters(filters, width_coefficient, depth_divisor):
    """Round number of filters based on width multiplier."""

    filters *= width_coefficient
    new_filters = int(filters + depth_divisor / 2) // depth_divisor * depth_divisor
    new_filters = max(depth_divisor, new_filters)
    # Make sure that round down does not go down by more than 10%.
    if new_filters < 0.9 * filters:
        new_filters += depth_divisor
    return int(new_filters)


def round_repeats(repeats, depth_coefficient):
    """Round number of repeats based on depth multiplier."""

    return int(math.ceil(depth_coefficient * repeats))



def EfficientNet_1d(width_coefficient,
                 depth_coefficient,
                 default_resolution,
                 dropout_rate=0.2,
                 drop_connect_rate=0.2,
                 depth_divisor=8,
                 blocks_args=DEFAULT_BLOCKS_ARGS,
                 model_name='efficientnet',
                 include_top=True,
                 input_tensor=None,
                 input_shape=None,
                 pooling=None,
                 output_shape=20,
                 activation = tf.nn.relu6):
    """Instantiates the EfficientNet architecture using given scaling coefficients.
    Optionally loads weights pre-trained on ImageNet.
    Note that the data format convention used by the model is
    the one specified in your Keras config at `~/.keras/keras.json`.
    # Arguments
        width_coefficient: float, scaling coefficient for network width.
        depth_coefficient: float, scaling coefficient for network depth.
        default_resolution: int, default input image size.
        dropout_rate: float, dropout rate before final classifier layer.
        drop_connect_rate: float, dropout rate at skip connections.
        depth_divisor: int.
        blocks_args: A list of BlockArgs to construct block modules.
        model_name: string, model name.
        include_top: whether to include the fully-connected
            layer at the top of the network.
        weights: one of `None` (random initialization),
              'imagenet' (pre-training on ImageNet),
              or the path to the weights file to be loaded.
        input_tensor: optional Keras tensor
            (i.e. output of `layers.Input()`)
            to use as image input for the model.
        input_shape: optional shape tuple, only to be specified
            if `include_top` is False.
            It should have exactly 3 inputs channels.
        pooling: optional pooling mode for feature extraction
            when `include_top` is `False`.
            - `None` means that the output of the model will be
                the 4D tensor output of the
                last convolutional layer.
            - `avg` means that global average pooling
                will be applied to the output of the
                last convolutional layer, and thus
                the output of the model will be a 2D tensor.
            - `max` means that global max pooling will
                be applied.
        classes: optional number of classes to classify images
            into, only to be specified if `include_top` is True, and
            if no `weights` argument is specified.
    # Returns
        A Keras model instance.
    # Raises
        ValueError: in case of invalid argument for `weights`,
            or invalid input shape.
    """

    # Determine proper input shape
    X_input = layers.Input(shape = (input_shape, 1))

    bn_axis = 2

    # Build stem
    x = X_input
    x = layers.Conv1D(round_filters(32, width_coefficient, depth_divisor), 3,
                      strides=2,
                      padding='same',
                      use_bias=False,
                      kernel_initializer=CONV_KERNEL_INITIALIZER,
                      name='stem_conv')(x)
    x = layers.BatchNormalization(axis=bn_axis, name='stem_bn')(x)
    x = layers.Activation(activation, name='stem_activation')(x)

    # Build blocks
    num_blocks_total = sum(round_repeats(block_args.num_repeat,
                                         depth_coefficient) for block_args in blocks_args)
    block_num = 0
    for idx, block_args in enumerate(blocks_args):
        assert block_args.num_repeat > 0
        # Update block input and output filters based on depth multiplier.
        block_args = block_args._replace(
            input_filters=round_filters(block_args.input_filters,
                                        width_coefficient, depth_divisor),
            output_filters=round_filters(block_args.output_filters,
                                         width_coefficient, depth_divisor),
            num_repeat=round_repeats(block_args.num_repeat, depth_coefficient))

        # The first block needs to take care of stride and filter size increase.
        drop_rate = drop_connect_rate * float(block_num) / num_blocks_total
        x = mb_conv_block_1d(x, block_args,
                          activation=activation,
                          drop_rate=drop_rate,
                          prefix='block{}a_'.format(idx + 1))
        block_num += 1
        if block_args.num_repeat > 1:
            # pylint: disable=protected-access
            block_args = block_args._replace(
                input_filters=block_args.output_filters, strides=1)
            # pylint: enable=protected-access
            for bidx in xrange(block_args.num_repeat - 1):
                drop_rate = drop_connect_rate * float(block_num) / num_blocks_total
                block_prefix = 'block{}{}_'.format(
                    idx + 1,
                    string.ascii_lowercase[bidx + 1]
                )
                x = mb_conv_block_1d(x, block_args,
                                  activation=activation,
                                  drop_rate=drop_rate,
                                  prefix=block_prefix)
                block_num += 1

    # Build top
    x = layers.Conv1D(round_filters(1280, width_coefficient, depth_divisor), 1,
                      padding='same',
                      use_bias=False,
                      kernel_initializer=CONV_KERNEL_INITIALIZER,
                      name='top_conv')(x)
    x = layers.BatchNormalization(axis=bn_axis, name='top_bn')(x)
    x = layers.Activation(activation, name='top_activation')(x)
    
    if include_top:
        x = layers.GlobalAveragePooling1D(name='avg_pool')(x)
        if dropout_rate and dropout_rate > 0:
            x = layers.Dropout(dropout_rate, name='top_dropout')(x)
        x = Dense(output_shape, name='output')(x)
    else:
        if pooling == 'avg':
            x = layers.GlobalAveragePooling1D(name='avg_pool')(x)
        elif pooling == 'max':
            x = layers.GlobalMaxPooling1D(name='max_pool')(x)

    # Create model.
    model = Model(inputs= X_input, outputs = x, name=model_name)
    
    return model



def EfficientNet_1dB0(input_shape,pooling='max', include_top=True,output_shape=20):
    return EfficientNet_1d(
        1.0, 1.0, 224, 0.2,
        model_name='efficientnet-b0',
        include_top=include_top,
        input_shape=input_shape,
        pooling=pooling, output_shape=output_shape
    )


def EfficientNet_1dB1(input_shape,pooling='max', include_top=True,output_shape=20):
    return EfficientNet_1d(
        1.0, 1.1, 240, 0.2,
        model_name='efficientnet-b1',
        include_top=include_top,
        input_shape=input_shape,
        pooling=pooling, output_shape=output_shape
    )


def EfficientNet_1dB2(input_shape,pooling='max', include_top=True,output_shape=20):
    return EfficientNet_1d(
        1.1, 1.2, 260, 0.3,
        model_name='efficientnet-b2',
        include_top=include_top,
        input_shape=input_shape,
        pooling=pooling, output_shape=output_shape
    )


def EfficientNet_1dB3(input_shape,pooling='max', include_top=True,output_shape=20):
    return EfficientNet_1d(
        1.2, 1.4, 300, 0.3,
        model_name='efficientnet-b3',
        include_top=include_top,
        input_shape=input_shape,
        pooling=pooling, output_shape=output_shape
    )


def EfficientNet_1dB4(input_shape,pooling='max', include_top=True,output_shape=20):
    return EfficientNet_1d(
        1.4, 1.8, 380, 0.4,
        model_name='efficientnet-b4',
        include_top=include_top,
        input_shape=input_shape,
        pooling=pooling, output_shape=output_shape
    )


def EfficientNet_1dB5(input_shape,pooling='max', include_top=True,output_shape=20):
    return EfficientNet_1d(
        1.6, 2.2, 456, 0.4,
        model_name='efficientnet-b5',
        include_top=include_top,
        input_shape=input_shape,
        pooling=pooling, output_shape=output_shape
    )


def EfficientNet_1dB6(input_shape,pooling='max', include_top=True,output_shape=20):
    return EfficientNet_1d(
        1.8, 2.6, 528, 0.5,
        model_name='efficientnet-b6',
        include_top=include_top,
        input_shape=input_shape,
        pooling=pooling, output_shape=output_shape
    )


def EfficientNet_1dB7(input_shape,pooling='max', include_top=True,output_shape=20):
    return EfficientNet_1d(
        2.0, 3.1, 600, 0.5,
        model_name='efficientnet-b7',
        include_top=include_top,
        input_shape=input_shape,
        pooling=pooling, output_shape=output_shape
    )


def EfficientNet_1dL2(input_shape,pooling='max', include_top=True,output_shape=20):
    return EfficientNet_1d(
        4.3, 5.3, 800, 0.5,
        model_name='efficientnet-l2',
        include_top=include_top,
        input_shape=input_shape,
        pooling=pooling, output_shape=output_shape
    )
