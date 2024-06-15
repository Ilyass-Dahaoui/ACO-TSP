import numpy as np
import matplotlib.pyplot as plt
import time
from tkinter import *
from tkinter import ttk
import tkinter as tk
from PIL import ImageTk,Image
import random as rd
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
plt.close()

pheromone_quantity = 100.0
alpha = 1.0
beta = 0.1
rho = 0.5 
number_of_ants = 50
number_of_steps = 100

class Ant: 
    def __init__(self, pheromone_quantity, alpha, beta, rho, graph, start):
        self.pheromone_quantity = pheromone_quantity
        self.alpha = alpha  
        self.beta = beta  
        self.rho = rho  
        self.carry_food = 0 
        self.distance_travelled = 0.0  
        self.node_road = [graph.nest]  
        self.edge_road = []  
        self.position = graph.nest 
        self.start = start
        self.graph = graph
        return
    
    def get_possible_edges(self):
        self.possible_edges = []
        for edge in self.graph.edges:
            if edge.start_node == self.position:
                self.possible_edges.append(edge)
        return self.possible_edges
    
    def select_edge(self):
        if self.start == 1:
            tau = np.array([edge.pheromone_level for edge in self.possible_edges])
            eta = np.array([1 / edge.length for edge in self.possible_edges])
            if np.sum(tau) == 0:
                self.best_edge = self.possible_edges[rd.randint(0, len(self.possible_edges) - 1)]
            else:
                p_edge = ((tau ** self.alpha) * (eta ** self.beta)) / np.sum((tau ** self.alpha) * (eta ** self.beta))
                index = np.argmax(p_edge)
                self.best_edge = self.possible_edges[index]
        else: 
            self.best_edge = self.possible_edges[rd.randint(0, len(self.possible_edges) - 1)]                                       
        return self.best_edge
    
    def update_position(self):
        self.node_road.append(self.best_edge.end_node)
        self.edge_road.append(self.best_edge)
        self.position = self.best_edge.end_node
        self.distance_travelled += self.best_edge.length
        return
    
    def pick_up_food(self):
        if self.position == self.graph.food:
            self.carry_food = 1
            self.start = 1
        return self.carry_food
    
    def update_pheromone(self):
        L = []
        for edge in self.graph.edges:
            if edge in self.edge_road:
                edge.pheromone_level += self.pheromone_quantity / self.distance_travelled
            edge.pheromone_level += (1 - self.rho) * edge.pheromone_level
            L.append(edge.pheromone_level)
        self.graph.max_phero = max(L)
        return
    
    def return_home(self):
        if self.carry_food == 1:
            self.position = self.graph.nest
            self.node_road = [self.graph.nest]
            self.edge_road = []
            self.carry_food = 0
            self.distance_travelled = 0.0
        return
    
    def plot_road(self):
        for i in range(len(self.edge_road)):
            plt.plot([self.node_road[i].x, self.node_road[i + 1].x],
                     [self.node_road[i].y, self.node_road[i + 1].y], color='green', linewidth=2)
        return
    
    def ant_mainloop(self):
        while self.carry_food == 0:
            self.get_possible_edges()
            self.select_edge()
            self.update_position()
            self.pick_up_food()
        return
    
    def ant_step(self):
        if self.carry_food == 1:
            self.update_pheromone()
            self.return_home()
        else:
            self.get_possible_edges()
            self.select_edge()
            self.update_position()
            self.pick_up_food()
        return

class Colony:
    def __init__(self, number_of_ants, number_of_steps, pheromone_quantity, alpha, beta, rho, graph):
        self.number_of_ants = number_of_ants
        self.number_of_steps = number_of_steps
        self.pheromone_quantity = pheromone_quantity
        self.alpha = alpha 
        self.beta = beta 
        self.rho = rho 
        self.graph = graph
    
    def create_colony(self):
        self.ants = [Ant(self.pheromone_quantity, self.alpha, self.beta, self.rho, self.graph, 0) for _ in range(self.number_of_ants)]
        return self.ants
    
    def roam_colony(self):
        for ant in self.ants:
            ant.ant_mainloop()
        return
    
    def colony_step(self):
        self.graph.pheromone_plot()
        self.graph.plot_nodes()
        for ant in self.ants:
            ant.ant_step()
            self.graph.plot_ant(ant)
        self.graph.plot_food()
        self.graph.plot_nest()
        self.graph.graph_window()
        return
    
    def colony_round_trip(self):
        self.roam_colony()
        for ant in self.ants:
            ant.update_pheromone()
            ant.return_home()
        self.graph.pheromone_plot()
        self.graph.plot_nodes()
        self.graph.plot_food()
        self.graph.plot_nest()
        self.graph.graph_window()
        return
    
    def colony_mainloop(self):
        self.create_colony()
        self.PCC = []
        for _ in range(self.number_of_steps):
            self.roam_colony()
            self.PCC.append(self.ants[0].edge_road)
            for ant in self.ants:
                ant.update_pheromone()
                ant.return_home()
            self.graph.pheromone_plot()
            self.graph.plot_nodes()
            self.graph.plot_food()
            self.graph.plot_nest()
            self.graph.graph_fenetre()
            plt.pause(0.5)
        return
    
class Graph: 
    def __init__(self, graph, max_phero, intrfc):
        self.height = 10
        self.width = 5
        self.coord_nodes = self.normalize_coord(graph[0])
        self.id_nodes = graph[1]
        self.edges_id = graph[2]
        self.id_nest = graph[3]
        self.id_food = graph[4]
        self.total_number_of_nodes = len(self.id_nodes)
        self.fig = plt.figure(figsize=(self.height, self.width))
        self.max_phero = max_phero
        self.intrfc = intrfc
    
    def normalize_coord(self, coord_nodes):
        res = []
        for coords in coord_nodes:
            coord1, coord2 = coords
            res.append((coord1 / 100 - 4, -coord2 / 100 - 4))
        return res  
    
    def make_nodes(self):
        self.nodes = [Node(self.coord_nodes[i]) for i in range(self.total_number_of_nodes)]
        return self.nodes
    
    def plot_nodes(self):
        for i in range(len(self.nodes)):
            plt.scatter(self.nodes[i].x, self.nodes[i].y, s=80, marker='o', color='black')
        return
    
    def plot_ant(self, L):
        plt.scatter(L.position.x + (-1) ** (rd.randint(0, 1)) * rd.random() / 5, 
                    L.position.y + (-1) ** (rd.randint(0, 1)) * rd.random() / 5, s=10, marker='+', color='grey')
        return
    
    def make_edges(self): 
        self.edges = []
        for edge_id in self.edges_id:
            id1, id2 = edge_id
            index1 = self.id_nodes.index(id1)
            index2 = self.id_nodes.index(id2)
            self.edges.append(Edge(self.nodes[index1], self.nodes[index2]))
        return self.edges
    
    def plot_edges(self):
        for i in range(len(self.edges)):
            plt.plot([self.edges[i].start_node.x, self.edges[i].end_node.x],
                     [self.edges[i].start_node.y, self.edges[i].end_node.y], color='black', linewidth=0.1)
        return
    
    def choose_nest(self):
        self.nest = self.nodes[self.id_nodes.index(self.id_nest)]
        return self.nest
    
    def plot_nest(self):
        plt.scatter(self.nest.x, self.nest.y, color='green', label='nest')
        return
    
    def choose_food(self):
        self.food = self.nodes[self.id_nodes.index(self.id_food)]
        return self.food
    
    def plot_food(self):
        plt.scatter(self.food.x, self.food.y, color='yellow', label='Food')
        return
    
    def pheromone_plot(self):
        self.fig = plt.figure(figsize=(self.height, self.width))
        for i in range(len(self.edges)):
            if self.edges[i].pheromone_level == self.max_phero:
                plt.plot([self.edges[i].start_node.x, self.edges[i].end_node.x],
                         [self.edges[i].start_node.y, self.edges[i].end_node.y],
                         color='red', linewidth=0.1 + 4 * self.edges[i].pheromone_level / self.max_phero)
            else:
                plt.plot([self.edges[i].start_node.x, self.edges[i].end_node.x],
                         [self.edges[i].start_node.y, self.edges[i].end_node.y],
                         color='red', linewidth=0.1 + 2 * self.edges[i].pheromone_level / self.max_phero)
        return
               
    def graph_window(self):
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.intrfc) 
        self.canvas.get_tk_widget().grid(row=7, column=2)
        return
        
    def graph_mainloop(self):
        self.make_nodes()
        self.make_edges()
        self.plot_nodes()
        self.plot_edges()
        self.choose_nest()
        self.choose_food()
        self.plot_food()
        self.plot_nest()
        self.graph_window()
        return


class GraphicalInterface(Tk):
    def __init__(self, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)
        self.iconbitmap(default='question')
        self.geometry("600x400")
        self.title("Home window")
        container = Frame(self)
        container.pack(fill="both", expand=True)
        image = Image.open("nid.jpeg").resize((600, 400), Image.LANCZOS)
        self.photo = ImageTk.PhotoImage(image)
        background_label = Label(container, image=self.photo)
        background_label.pack(fill="both", expand=True)
        plot_button = Button(container, text="Ant Colony Optimization", command=self.destroy,
                             font=("Arial", 12), fg="white", bg="green", padx=10, pady=5, relief=RAISED)
        plot_button.place(relx=0.5, rely=0.5, anchor="center")
        return None
    
class DisplayCanvas(Canvas):
    def __init__(self, parent, w=500, h=400, _bg='white'):  
        self.__w = w
        self.__h = h
        self.__node_list = [] 
        self.__lines = [] 

        self.__parent_window = parent
        Canvas.__init__(self, parent, width=w, height=h, bg=_bg, relief=RAISED, bd=5)

    def get_dims(self):
        return (self.__w, self.__h)

    def not_used_redraw_graph(self):
        for node in self.__node_list:
            node.move()
        self.after(60, self.draw_graph)  

    def create_node(self, x_center, y_center, radius, color, fill_color="white"):
        node = Ball(self, x_center, y_center, radius, color)
        self.pack()
        return node

    def left_click(self, event):
        self.__parent_window.place_node(event.x, event.y)

    def display_path(self, path_list, node_list, node_id_list):
        for path in path_list:
            first_node = node_list[node_id_list.index(path[0])]
            second_node = node_list[node_id_list.index(path[1])]
            self.__lines.append(self.create_line(first_node, second_node, fill="red", width=1))

    def erase_paths(self):
        for line in self.__lines:
            self.delete(line)

    def not_used_set_last_node_coordinates(self, x_center, y_center):
        self.__last_node = (x_center, y_center)

    def place_node_on_canvas(self, x_center, y_center, color='black', fill_color="black"):
        w, h = self.get_dims()
        radius = 5
        node = self.create_node(x_center, y_center, radius, color, fill_color)
        self.update()
        self.__parent_window.set_last_node_coordinates(x_center, y_center)
        return node.get_node_id()
    
class MainWindow(Tk):
    def __init__(self):
        Tk.__init__(self)
        self.title("Please build your graph :)")
        self.__display_area = DisplayCanvas(self)
        self.__display_area.pack()
        self.current_path = False
        self.arrival = 0
        self.departure_id = -1
        self.arrival_id = -1
        self.text1 = StringVar()
        self.text1.set("Add paths")
        self.text2 = StringVar()
        self.text2.set("Arrival")
        self.__undo_button = Button(self, text='Undo', command=self.undo_last_node).pack(side=LEFT, padx=5, pady=5)
        self.__clear_button = Button(self, text='Clear', command=self.clear).pack(side=LEFT, padx=5, pady=5)
        self.__path_button = Button(self, textvariable=self.text1, command=self.toggle_current_path).pack(side=RIGHT, padx=5, pady=5)
        self.__arrival_button = Button(self, textvariable=self.text2, command=self.set_arrival_and_departure).pack(side=RIGHT, padx=5, pady=5)
        self.__save_button = Button(self, text='Save', command=self.save).pack(side=RIGHT, padx=5, pady=5)
        self.__quit_button = Button(self, text='Quit', command=self.destroy).pack(side=LEFT, padx=5, pady=5)
        self.__display_area.bind('<Button-1>', self.__display_area.left_click)
        self.__created_object_id_list = []
        self.__center_coordinates_list = []
        self.__current_path_construction = []
        self.__path_list = []

    def toggle_current_path(self):
        self.current_path = not self.current_path
        if self.current_path:
            self.__display_area.bind('<Button-1>', self.build_path)
            self.text1.set("Add nodes")
        else:
            self.__display_area.bind('<Button-1>', self.__display_area.left_click)
            self.text1.set("Add paths")
        return None

    def set_arrival_and_departure(self):
        self.arrival += 1
        if self.arrival == 3:
            self.arrival = 0
        if self.arrival == 0:
            self.__display_area.bind('<Button-1>', self.__display_area.left_click)
            self.text2.set("Arrival")
            print("Nodes")
        if self.arrival == 1:
            self.__display_area.bind('<Button-1>', self.add_arrival)
            self.text2.set("Departure")
            print("Arrival")
        if self.arrival == 2:
            self.__display_area.bind('<Button-1>', self.add_departure)
            self.text2.set("Nodes")
            print("Departure")
        return None

    def add_arrival(self, event):
        in_node = False
        for k in range(len(self.__center_coordinates_list)):
            x_center, y_center = self.__center_coordinates_list[k]
            if (event.x - x_center)**2 + (event.y - y_center)**2 < 25:
                node_id = self.__created_object_id_list[k]
                in_node = True
        if in_node:
            if self.arrival_id != -1:
                self.__display_area.itemconfig(self.arrival_id, fill="white")
            self.__display_area.itemconfig(node_id, fill="yellow")
            self.arrival_id = node_id

    def add_departure(self, event):
        in_node = False
        for k in range(len(self.__center_coordinates_list)):
            x_center, y_center = self.__center_coordinates_list[k]
            if (event.x - x_center)**2 + (event.y - y_center)**2 < 25:
                node_id = self.__created_object_id_list[k]
                in_node = True
        if in_node:
            if self.departure_id != -1:
                self.__display_area.itemconfig(self.departure_id, fill="white")
            self.__display_area.itemconfig(node_id, fill="green")
            self.departure_id = node_id

    def add_node_to_list(self, node) :
        self.__created_object_id_list.append(node)

    def place_node(self, x, y):
        node = self.__display_area.place_node_on_canvas(x, y)
        self.add_node_to_list(node)
        self.set_last_node_coordinates(x, y)
        self.__center_coordinates_list.append((x, y))

    def build_path(self, event):
        in_node = False
        for k in range(len(self.__center_coordinates_list)):
            x_center, y_center = self.__center_coordinates_list[k]
            if (event.x - x_center)**2 + (event.y - y_center)**2 < 25:
                node_id = self.__created_object_id_list[k]
                in_node = True
        if in_node:
            if self.__current_path_construction != []:
                if self.__current_path_construction[0] != node_id:
                    self.__current_path_construction.append(node_id)
            else:
                self.__current_path_construction.append(node_id)
                print("First Node")
            if len(self.__current_path_construction) == 2:
                print("Second Node")
                first_node = self.__current_path_construction[0]
                if (first_node, node_id) in self.__path_list or (node_id, first_node) in self.__path_list:
                    print("Already existing path")
                    self.__current_path_construction = []
                else:
                    self.__path_list.append((first_node, node_id))
                    self.__current_path_construction = []
                    self.__display_area.display_path(self.__path_list, self.__center_coordinates_list, self.__created_object_id_list)

    def set_last_node_coordinates(self, x_center, y_center):
        self.__last_node = (x_center, y_center)

    def get_last_node(self):
        return self.__last_node

    def undo_last_node(self):
        if len(self.__created_object_id_list) == 0:
            print("No Nodes to remove")
            return
        node_num = self.__created_object_id_list[-1]
        x_center, y_center = self.get_last_node()
        paths_to_remove = []
        for path in self.__path_list:
            x, y = path
            if x == node_num or y == node_num:
                paths_to_remove.append(path)
        for path in paths_to_remove:
            self.__path_list.remove(path)

        last_node = self.__created_object_id_list.pop()
        self.__display_area.delete(last_node)
        self.__display_area.update()
        x_last_node, y_last_node = self.__center_coordinates_list.pop()
        self.set_last_node_coordinates(x_last_node, y_last_node)
        self.__display_area.erase_paths()
        self.__display_area.display_path(self.__path_list, self.__center_coordinates_list, self.__created_object_id_list)

    def add_node(self):
        x_center,  y_center = self.not_used_generate_random_point()
        radius = 5
        color = 'black'
        self.__display_area.create_node(x_center, y_center, radius, color)
        self.__display_area.update()
        self.__last_node = (x_center, y_center)

    def clear(self):
        self.__display_area.delete(ALL)
        self.__created_object_id_list.clear()
        self.__center_coordinates_list.clear()
        self.__path_list.clear()

    def save(self):
        if self.departure_id == -1:
            print("Departure not defined")
            return
        if self.arrival_id == -1:
            print("Arrival not defined")
            return
        print("Saved graph")
        self.graph = self.get_graph()
        self.destroy()

    def get_graph(self):
        return (self.__center_coordinates_list, 
                self.__created_object_id_list,
                self.__path_list,
                self.departure_id,
                self.arrival_id)

class Ball:
    def __init__(self, canvas, cx, cy, radius, color, fill_color="black"):
        self.__cx, self.__cy = cx, cy
        self.__radius = radius
        self.__color = color
        self.__canvas = canvas  
        self.__canvas_id = self.__canvas.create_oval(cx - radius, cy - radius, cx + radius, cy + radius, outline=color, fill=fill_color)

    def get_node_id(self):
        return self.__canvas_id

    def move(self):
        self.__canvas.move(self.__canvas_id, 0, 0)
        
class App(Tk):
    def __init__(self, graph, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)
        container = Frame(self)
        self.iconbitmap(default='warning')
        self.title("Main Program")
        self.geometry("600x400")
        self.init_graph = graph
        
        self.apha= alpha
        self.beta = beta
        self.rho = rho
        self.number_of_ants = number_of_ants

        # Titre
        Label(self, text="ACO Parameters", font=("Helvetica", 10)).grid(row=0, column=0, columnspan=2, pady=10)

        # Configuration des paramètres
        Label(self, text="Number of ants:").grid(row=1, column=0, padx=10, pady=5, sticky="e")
        self.number_of_ants_entry = Entry(self, width=10)
        self.number_of_ants_entry.grid(row=1, column=1, padx=10, pady=5)
        self.number_of_ants_entry.insert(0,str(self.number_of_ants))


        Label(self, text="Alpha:").grid(row=2, column=0, padx=10, pady=5, sticky="e")
        self.alpha_entry = Entry(self, width=10)
        self.alpha_entry.grid(row=2, column=1, padx=10, pady=5)
        self.alpha_entry.insert(0,str(self.apha))


        Label(self, text="Beta:").grid(row=3, column=0, padx=10, pady=5, sticky="e")
        self.beta_entry = Entry(self, width=10)
        self.beta_entry.grid(row=3, column=1, padx=10, pady=5)
        self.beta_entry.insert(0,str(self.beta))

        Label(self, text="Rho:").grid(row=4, column=0, padx=10, pady=5, sticky="e")
        self.rho_entry = Entry(self, width=10)
        self.rho_entry.grid(row=4, column=1, padx=10, pady=5)
        self.rho_entry.insert(0,str(self.rho))

        # Bouton de démarrage
        self.update_button = Button(self, command=self.start, text="Start", bg="green", fg="white", width=10)
        self.update_button.grid(row=5, column=0, columnspan=2, pady=20)

        # Ajustement des dimensions des colonnes et lignes
        self.grid_columnconfigure(1, weight=1) 
        self.grid_rowconfigure(5, weight=1)

   

    def start(self):
        self.alpha = float(self.alpha_entry.get())
        self.beta = float(self.beta_entry.get())
        self.rho = float(self.rho_entry.get())
        self.number_of_ants = int(self.number_of_ants_entry.get())
        self.graph = Graph(self.init_graph, 0.9, app)
        self.graph.graph_mainloop()
        self.colony = Colony(self.number_of_ants, number_of_steps, pheromone_quantity, self.alpha, self.beta, self.rho, self.graph)
        self.colony.create_colony()
        self.make_round_trip_button = Button(self, command=self.colony.colony_round_trip, text="Make a round trip") 
        self.make_round_trip_button.grid(row=50, column=2)
        self.take_step_button = Button(self, command=self.colony.colony_step, text="Take a step") 
        self.take_step_button.grid(row=51, column=2)
        
class Node:     
    def __init__(self, coord):
        self.x = coord[0]
        self.y = coord[1]

class Edge:     
    def __init__(self, start_node, end_node):
        self.length = np.sqrt((start_node.x - end_node.x) ** 2 + (start_node.y - end_node.y) ** 2)
        self.start_node = start_node
        self.end_node = end_node
        self.pheromone_level = 0.1


        
if __name__ == "__main__":
    # Creating and running the InterfaceGraphique page
    page = GraphicalInterface()
    page.mainloop()

    # Creating and running the MainWindow
    fen = MainWindow()
    fen.mainloop()

    # Getting the graph from the MainWindow
    graphe = fen.get_graph()

    # Checking if a graph was retrieved successfully
    if graphe != ([], [], [], -1, -1):
        # Creating and running the App with the retrieved graph
        app = App(graphe)
        app.mainloop()
        
