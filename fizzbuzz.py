output = ''
for i in range(1, 101):
    #print(i)
    if i % 3 == 0:
        output = f'fizz{output}'
    if i % 5 == 0:
        output = f'{output}buzz'
        print(output, i)
        output = ''
    if len(output) > 1:
        print(output, i)
        output = ''