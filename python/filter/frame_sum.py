import tensorflow as tf
import numpy as np

def frameSum(arr):
    
    x_data = tf.placeholder(tf.uint16, shape=(1742000))

    X = tf.Variable(x_data)

    frameSum = tf.reduce_sum(x_data)

    with tf.Session() as session:
        tf.initialize_all_variables()
        res = session.run(frameSum, feed_dict={x_data: arr})

    return res
