import matplotlib.pyplot as plt

# Data based on the 33 ligands
classes = [
    'Flavonoids', 
    'Phenolic Acids', 
    'Triterpenoids/Sterols', 
    'Sesquiterpenes', 
    'N-containing Compounds', 
    'Fatty Acids/Amides', 
    'Others'
]
counts = [7, 5, 7, 5, 2, 2, 5] # Adjusted counts based on the list above

# Colors for each slice
colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99', '#c2c2f0', '#ffb3e6', '#c4e17f']

# Create Pie Chart
plt.figure(figsize=(10, 7))
plt.pie(counts, labels=classes, autopct='%1.1f%%', startangle=140, colors=colors, pctdistance=0.85)

# Draw a circle at the center to make it a donut chart (optional, looks more modern)
centre_circle = plt.Circle((0,0),0.70,fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)

plt.title('Chemical Class Distribution of 33 Dicoma anomala Ligands', fontsize=14)
plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.tight_layout()
plt.savefig('chemical_class_distribution.png', dpi=300)
plt.show()