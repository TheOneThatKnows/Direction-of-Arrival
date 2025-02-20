import tensorflow as tf
from tensorflow.keras import layers, Model
import numpy as np

class ComplexDense(layers.Layer):
    def __init__(self, units):
        super(ComplexDense, self).__init__()
        self.units = units
        
    def build(self, input_shape):
        # input_shape will be (batch_size, input_dim, 2)
        input_dim = input_shape[-2]
        
        self.w_real = self.add_weight(
            name='w_real',
            shape=(input_dim, self.units),
            initializer='glorot_uniform',
            trainable=True)
        
        self.w_imag = self.add_weight(
            name='w_imag',
            shape=(input_dim, self.units),
            initializer='glorot_uniform',
            trainable=True)
        
        self.b_real = self.add_weight(
            name='b_real',
            shape=(self.units,),
            initializer='zeros',
            trainable=True)
        
        self.b_imag = self.add_weight(
            name='b_imag',
            shape=(self.units,),
            initializer='zeros',
            trainable=True)
        
    def call(self, inputs):
        # inputs shape: (batch_size, input_dim, 2)
        real_part = inputs[..., 0]  # Shape: (batch_size, input_dim)
        imag_part = inputs[..., 1]  # Shape: (batch_size, input_dim)
        
        # Perform matrix multiplication
        output_real = tf.matmul(real_part, self.w_real) - tf.matmul(imag_part, self.w_imag) + self.b_real
        output_imag = tf.matmul(real_part, self.w_imag) + tf.matmul(imag_part, self.w_real) + self.b_imag
        
        # Stack real and imaginary parts
        return tf.stack([output_real, output_imag], axis=-1)

def create_model(input_shape, output_shape):
    # Input layer
    inputs = layers.Input(shape=input_shape)
    
    # Initial normalization
    x = layers.BatchNormalization()(inputs)
    
    # First block
    x1 = ComplexDense(256)(x)
    x1 = layers.LeakyReLU(alpha=0.2)(x1)
    x1 = layers.Dropout(0.2)(x1)
    
    # Skip connection block 1
    x2 = ComplexDense(256)(x1)
    x2 = layers.LeakyReLU(alpha=0.2)(x2)
    x2 = layers.Dropout(0.2)(x2)
    x2 = layers.Add()([x1, x2])  # Skip connection
    
    # Skip connection block 2
    x3 = ComplexDense(128)(x2)
    x3 = layers.LeakyReLU(alpha=0.2)(x3)
    x3 = layers.Dropout(0.2)(x3)
    
    # Skip connection block 3
    x4 = ComplexDense(128)(x3)
    x4 = layers.LeakyReLU(alpha=0.2)(x4)
    x4 = layers.Add()([x3, x4])  # Skip connection
    
    # Final blocks
    x5 = ComplexDense(64)(x4)
    x5 = layers.LeakyReLU(alpha=0.2)(x5)
    
    # Output layer
    outputs = ComplexDense(output_shape)(x5)
    
    # Create model
    model = Model(inputs=inputs, outputs=outputs)
    
    return model