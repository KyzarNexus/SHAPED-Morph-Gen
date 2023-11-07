"""
GUI_TEST.py
Just tinkering around with the idea of a basic GUI to manipulate and test data and more succinctly display certain graphs. 
"""
#%%
import tkinter
from tkinter import messagebox

top = tkinter.Tk()

def helloCallBack():
   messagebox.showinfo( "Hello Python", "Hello World")

B = tkinter.Button(top, text ="Hello", command = helloCallBack)

B.pack()
top.mainloop()
# %%




















