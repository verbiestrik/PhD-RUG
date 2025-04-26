import configparser

config = configparser.ConfigParser()

config.read('/Users/vildarozek/Files/Fork/Untitled/Spatially Adaptive Moment Models/Nonlinear-systems/Config-files/config.txt')
print(config.sections())

numerical_method_information = config['numerical_method_information']

ordersList = numerical_method_information['orders']
print(ordersList)
ordersList = [int(order) for order in ordersList.split(',')]

print(ordersList)

