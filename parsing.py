#!/usr/bin/env python
# coding: utf-8

# In[11]:


import os
os.chdir(r'C:\Users\amand\Downloads')

def parse(filename):
    with open(filename, 'r') as handle:
        #handle = 
        handle = handle.read().split('\n')
        key =[]
        value = []
        for line in handle:
            #print(line)
            if line.startswith('>>Basic Statistics'):
                a = str((line[2:]))
                a = a.split('\t')
                key.append(a[0])
                value.append(a[1])

            if line.startswith('Filename'):
                a = str(line)
                a = a.split('\t')
                key.append(a[0])
                value.append(a[1])
            if line.startswith('Total Sequences'):
                a = str(line)
                a = a.split('\t')
                key.append(a[0])
                value.append(a[1])
            if line.startswith('>>Sequence Duplication Levels'):
                a = str(line[2:])
                a = a.split('\t')
                key.append(a[0])
                value.append(a[1])
            if line.startswith('#Total Deduplicated Percentage'):
                a = str(line[1:])
                a = a.split('\t')
                key.append(a[0])
                b= round(int(float(a[1])))
                per = str(b) + "%"
                value.append(per)
            if line.startswith('>>Overrepresented sequences'):
                a = str(line[2:])
                a = a.split('\t')
                key.append(a[0])
                value.append(a[1])
            if line.startswith('>>Adapter Content'):
                a = str(line[2:])
                a = a.split('\t')
                key.append(a[0])
                value.append(a[1])

    #print(key)
    #print(value)

        my_dictionary = dict(zip(key, value))



        print(my_dictionary)
    
parse("fastqc_data.txt")




  
    

    



  


# In[ ]:




