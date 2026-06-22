import plotly.graph_objects as go

def init_plot():
    """ 
    Initialize the plot 
    """
    fig = go.Figure()
    fig.add_trace(
        go.Scatter3d(x=[0], y=[0], z=[0], mode = 'markers', name = "Black hole",
            marker = dict(
                color='black',
                size=10
                ))
            )
    return fig

def plot(fig, coords, label=None):
    """ 
    Plot some orbit/trajectory on specified plot 
    """
    x, y, z = coords
    show_label = True if label else False
    
    fig.add_trace(go.Scatter3d(
        x = x, y = y, z = z,
        mode = 'lines',
        name = label,
        showlegend = show_label
    ))
    
def show_plot(fig):
    """
    Show the plot and save it as HTML file
    """
    fig.update_layout(
        scene=dict(
            xaxis_title = 'x',
            yaxis_title = 'y',
            zaxis_title = 'z',
            aspectmode='data'
        ),
        title = "TDE Simulation in PN Approximation"
    )
    fig.show()
    fig.write_html('./graphics/tde_3d.html')