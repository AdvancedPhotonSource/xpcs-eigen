import tensorflow as tf
import numpy as np

def pixelSum(arr):

    x_data = tf.placeholder(tf.uint16, shape=(266, 1742000))

    X = tf.Variable(x_data)

    pixelSum = tf.reduce_sum(x_data, 0)

    with tf.Session() as session:
        tf.initialize_all_variables()
        res = session.run(pixelSum, feed_dict={x_data: arr})

    return res
