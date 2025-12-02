# trace generated using paraview version 5.13.3
# import paraview
# paraview.compatibility.major = 5
# paraview.compatibility.minor = 13

import glob
import os
import re

#### import the simple module from the paraview
from paraview.simple import *

#### disable automatic camera reset on 'Show'

sod_test_n = ["sod1", "sod2", "sod3", "sod4", "sod5"]

legend_loc = [
    ["TopRight", "TopLeft", "TopRight", "Custom"],  # 1
    ["Custom", "TopLeft", "Custom", "Custom"],  # 2
    ["TopLeft", "TopRight", "TopRight", "TopRight"],  # 3
    ["TopRight", "Custom", "TopLeft", "TopLeft"],  # 4
    ["TopLeft", "TopRight", "TopLeft", "TopLeft"]  # 5
]

line_style_dict = {
    "analytical": "2",
    "godunov": "1",
    "godunov-kolgan": "1",
    "godunov-kolgan-rodionov": "1"
}

label_dict = {
    "analytical": "Analytical",
    "godunov": "Godunov",
    "godunov-kolgan": "Godunov-Kolgan",
    "godunov-kolgan-rodionov": "Godunov-Kolgan-Rodionov"
}

# Функция для извлечения солвера и реконструкции из имени папки
def parse_folder_name(folder_name):
    if folder_name == "analytical":
        return "analytical", ""

    # Регулярное выражение для извлечения солвера и реконструкции
    match = re.match(r'^([a-z-]+)__R_([a-z0-9]+)__N_\d+__CFL_[\d.e-]+$', folder_name)
    if match:
        solver = match.group(1)
        reconstruction = match.group(2)
        return solver, reconstruction
    return None, None

# Функция для генерации уникального цвета
def get_color(solver, reconstruction, color_index):
    if solver == "analytical":
        return ["0", "0", "0"]  # Черный для analytical

    # Генерация различных цветов для разных реконструкций
    colors = [
        ["0", "0", "1"],       # Синий
        ["1", "0", "0"],       # Красный
        ["0", "1", "0"],       # Зеленый
        ["1", "1", "0"],       # Желтый
        ["1", "0", "1"],       # Пурпурный
        ["0", "1", "1"],       # Голубой
        ["0.5", "0", "0"],     # Темно-красный
        ["0", "0.5", "0"],     # Темно-зеленый
        ["0", "0", "0.5"],     # Темно-синий
        ["0.5", "0.5", "0"],   # Оливковый
        ["1", "0.5", "0"],     # Оранжевый
        ["0.5", "0", "0.5"],   # Фиолетовый
        ["0", "0.5", "0.5"],   # Темно-бирюзовый
        ["0.75", "0.25", "0"], # Коричневатый
        ["0.25", "0.75", "0"], # Салатовый
        ["0", "0.25", "0.75"], # Сине-голубой
        ["0.75", "0", "0.25"], # Розово-красный
        ["0.25", "0", "0.75"], # Темно-фиолетовый
        ["0.75", "0.75", "0"], # Золотистый
        ["0", "0.75", "0.75"], # Бирюзовый
        ["0.75", "0", "0.75"], # Ярко-фиолетовый
        ["0.3", "0.7", "0.3"], # Светло-зеленый
        ["0.7", "0.3", "0.3"], # Кирпичный
        ["0.3", "0.3", "0.7"], # Светло-синий
        ["0.8", "0.2", "0.5"], # Розовый
        ["0.2", "0.8", "0.5"], # Мятный
        ["0.5", "0.2", "0.8"], # Лавандовый
        ["0.8", "0.5", "0.2"], # Персиковый
        ["0.2", "0.5", "0.8"], # Небесный
        ["0.6", "0.4", "0.2"], # Карамельный
    ]

    return colors[color_index % len(colors)]

base_dir = "data"

# Основной цикл по тестам
for test_num, test_name in enumerate(sod_test_n):
    base_path = os.path.join(base_dir, test_name)

    # Проверяем существование папки теста
    if not os.path.exists(base_path):
        print(f"Warning: Test folder {base_path} not found")
        continue

    # Получаем все подпапки в папке теста
    all_folders = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))]

    if not all_folders:
        print(f"Warning: No subfolders found in {base_path}")
        continue

    print(f"Processing test {test_name}, found {len(all_folders)} folders")

    # Создаем новый layout для каждого теста
    layout = CreateLayout(name=f'Sod test {test_num + 1}')

    # Создаем 4 графика для каждого layout с правильным разделением
    lineChartViews = []

    # Создаем первый view
    lineChartView1 = CreateView('XYChartView')
    AssignViewToLayout(view=lineChartView1, layout=layout, hint=0)
    lineChartViews.append(lineChartView1)

    # Разделяем первую ячейку по горизонтали
    layout.SplitHorizontal(0, 0.5)

    # Создаем второй view в правой части
    lineChartView2 = CreateView('XYChartView')
    AssignViewToLayout(view=lineChartView2, layout=layout, hint=2)
    lineChartViews.append(lineChartView2)

    # Разделяем левую часть (ячейка 0) по вертикали
    layout.SplitVertical(1, 0.5)

    # Создаем третий view в левой нижней части
    lineChartView3 = CreateView('XYChartView')
    AssignViewToLayout(view=lineChartView3, layout=layout, hint=4)
    lineChartViews.append(lineChartView3)

    # Разделяем правую часть (ячейка 2) по вертикали
    layout.SplitVertical(2, 0.5)

    # Создаем четвертый view в правой нижней части
    lineChartView4 = CreateView('XYChartView')
    AssignViewToLayout(view=lineChartView4, layout=layout, hint=6)
    lineChartViews.append(lineChartView4)

    # Переименовываем views
    view_names = ['Density', 'Velocity', 'Pressure', 'Energy']
    for view, name in zip(lineChartViews, view_names):
        SetActiveView(view)
        RenameView(name, view)

    # Счетчик для цветов
    color_counter = 0

    # Обрабатываем каждую папку
    for folder_name in sorted(all_folders):
        solver, reconstruction = parse_folder_name(folder_name)

        if solver is None:
            print(f"Warning: Could not parse folder name: {folder_name}")
            continue

        folder_path = os.path.join(base_path, folder_name)

        # Ищем файлы .vtk
        if solver == "analytical":
            files = sorted(glob.glob(os.path.join(folder_path, 'step_*.vtk')))
        else:
            files = sorted(glob.glob(os.path.join(folder_path, f'{folder_name}__step_*.vtk')))

        if not files:
            print(f"Warning: No .vtk files found in {folder_path}")
            continue

        print(f"  Processing: {solver} ({reconstruction}) - {len(files)} files")

        # создаем новый 'Legacy VTK Reader'
        vtk = LegacyVTKReader(registrationName=f'VTK_test{test_num}_{solver}_{reconstruction}', FileNames=files)

        # get animation scene
        animationScene1 = GetAnimationScene()

        # update animation scene based on data timesteps
        animationScene1.UpdateAnimationUsingDataTimeSteps()

        # set active source
        SetActiveSource(vtk)

        animationScene1.GoToLast()

        # Получаем цвет для текущего набора данных
        color = get_color(solver, reconstruction, color_counter)
        color_counter += 1

        # Цикл по графикам (4 типа)
        plot_variables = ['density', 'velocity', 'pressure', 'conserved_energy']

        for i, (view, variable) in enumerate(zip(lineChartViews, plot_variables)):
            # создаем новый 'Plot Over Line'
            plotOverLine = PlotOverLine(
                registrationName=f'PlotOverLine_test{test_num}_{solver}_{reconstruction}_{variable}',
                Input=vtk)

            # Properties modified on plotOverLine
            plotOverLine.Point1 = [0.0005000000237487257, 0.49950000643730164, 0.0]
            plotOverLine.Point2 = [0.9994999766349792, 0.49950000643730164, 0.0]
            plotOverLine.Resolution = 1000

            # show data in view
            plotOverLineDisplay = Show(plotOverLine, view, 'XYChartRepresentation')

            # update the view to ensure updated data information
            view.Update()

            # Настройки отображения
            plotOverLineDisplay.SeriesOpacity = [
                'arc_length', '1',
                'conserved_energy', '1',
                'density', '1',
                'mass', '1',
                'momentum', '1',
                'pressure', '1',
                'specific_internal_energy', '1',
                'velocity', '1',
                'volume', '1',
                'vtkValidPointMask', '1',
                'Points_X', '1',
                'Points_Y', '1',
                'Points_Z', '1',
                'Points_Magnitude', '1'
            ]
            plotOverLineDisplay.SeriesPlotCorner = [
                'Points_Magnitude', '0',
                'Points_X', '0',
                'Points_Y', '0',
                'Points_Z', '0',
                'arc_length', '0',
                'conserved_energy', '0',
                'density', '0',
                'mass', '0',
                'momentum', '0',
                'pressure', '0',
                'specific_internal_energy', '0',
                'velocity', '0',
                'volume', '0',
                'vtkValidPointMask', '0'
            ]

            # Устанавливаем стиль линии
            line_style = line_style_dict.get(solver, "1")
            plotOverLineDisplay.SeriesLineStyle = [
                'Points_Magnitude', '1',
                'Points_X', '1',
                'Points_Y', '1',
                'Points_Z', '1',
                'arc_length', '1',
                'conserved_energy', line_style,
                'density', line_style,
                'mass', '1',
                'momentum', '1',
                'pressure', line_style,
                'specific_internal_energy', '1',
                'velocity', line_style,
                'volume', '1',
                'vtkValidPointMask', '1'
            ]

            plotOverLineDisplay.SeriesLineThickness = [
                'Points_Magnitude', '2',
                'Points_X', '2',
                'Points_Y', '2',
                'Points_Z', '2',
                'arc_length', '2',
                'conserved_energy', '2',
                'density', '2',
                'mass', '2',
                'momentum', '2',
                'pressure', '2',
                'specific_internal_energy', '2',
                'velocity', '2',
                'volume', '2',
                'vtkValidPointMask', '2'
            ]
            plotOverLineDisplay.SeriesMarkerStyle = [
                'Points_Magnitude', '0',
                'Points_X', '0',
                'Points_Y', '0',
                'Points_Z', '0',
                'arc_length', '0',
                'conserved_energy', '0',
                'density', '0',
                'mass', '0',
                'momentum', '0',
                'pressure', '0',
                'specific_internal_energy', '0',
                'velocity', '0',
                'volume', '0',
                'vtkValidPointMask', '0'
            ]
            plotOverLineDisplay.SeriesMarkerSize = [
                'Points_Magnitude', '4',
                'Points_X', '4',
                'Points_Y', '4',
                'Points_Z', '4',
                'arc_length', '4',
                'conserved_energy', '4',
                'density', '4',
                'mass', '4',
                'momentum', '4',
                'pressure', '4',
                'specific_internal_energy', '4',
                'velocity', '4',
                'volume', '4',
                'vtkValidPointMask', '4'
            ]

            # Устанавливаем видимость только для нужной переменной
            plotOverLineDisplay.SeriesVisibility = [variable]

            # Формируем метку
            if solver == "analytical":
                label = label_dict.get(solver, "Analytical")
            else:
                base_label = label_dict.get(solver, solver)
                label = f"{base_label} ({reconstruction})" if reconstruction else base_label

            # Устанавливаем метки
            plotOverLineDisplay.SeriesLabel = [
                'arc_length', 'arc_length',
                'conserved_energy', label,
                'density', label,
                'mass', 'mass',
                'momentum', 'momentum',
                'pressure', label,
                'specific_internal_energy', 'specific_internal_energy',
                'velocity', label,
                'volume', 'volume',
                'vtkValidPointMask', 'vtkValidPointMask',
                'Points_X', 'Points_X',
                'Points_Y', 'Points_Y',
                'Points_Z', 'Points_Z',
                'Points_Magnitude', 'Points_Magnitude']

            # Устанавливаем цвета
            plotOverLineDisplay.SeriesColor = [
                'arc_length', '0', '0', '0',
                'conserved_energy', color[0], color[1], color[2],
                'density', color[0], color[1], color[2],
                'mass', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155',
                'momentum', '0.6', '0.3100022888532845', '0.6399938963912413',
                'pressure', color[0], color[1], color[2],
                'specific_internal_energy', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867',
                'velocity', color[0], color[1], color[2],
                'volume', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845',
                'vtkValidPointMask', '0.220004577706569', '0.4899977111467155', '0.7199969481956207',
                'Points_X', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155',
                'Points_Y', '0.6', '0.3100022888532845', '0.6399938963912413',
                'Points_Z', '1', '0.5000076295109483', '0',
                'Points_Magnitude', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867'
            ]

            # Устанавливаем расположение легенды для каждого view
            legend_location = legend_loc[test_num][i]
            view.LegendLocation = legend_location
            if legend_location == 'Custom':
                view.LegendPosition = [view.ViewSize[0] // 2 - 50, 337]

    print(f"Completed processing test {test_name}")
